import logging
import re
import sqlite3
from pathlib import Path
from typing import Dict, Optional, Set

import click
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyopenms as oms

from quantmsutils.mzml.ms1_feature_finder import MS1FeatureDetector
from quantmsutils.openms import extract_scan_id
from quantmsutils.utils.constants import (
    ACQUISITION_DATETIME,
    BASE_PEAK_INTENSITY,
    CHARGE,
    EXPERIMENTAL_MASS_TO_CHARGE,
    INTENSITY,
    INTENSITY_ARRAY,
    MS_LEVEL,
    MZ_ARRAY,
    NUM_PEAKS,
    RETENTION_TIME,
    SCAN,
    SUMMED_PEAK_INTENSITY,
    PRECURSOR_RT,
    PRECURSOR_TOTAL_INTENSITY,
)

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


def create_ms_schema() -> pa.Schema:
    """Create and return the main spectra schema"""
    return pa.schema(
        [
            pa.field(SCAN, pa.string(), nullable=True),
            pa.field(MS_LEVEL, pa.int32(), nullable=True),
            pa.field(NUM_PEAKS, pa.int32(), nullable=True),
            pa.field(BASE_PEAK_INTENSITY, pa.float64(), nullable=True),
            pa.field(SUMMED_PEAK_INTENSITY, pa.float64(), nullable=True),
            pa.field(RETENTION_TIME, pa.float64(), nullable=True),
            pa.field(CHARGE, pa.int32(), nullable=True),
            pa.field(EXPERIMENTAL_MASS_TO_CHARGE, pa.float64(), nullable=True),
            pa.field(PRECURSOR_RT, pa.float64(), nullable=True),
            pa.field(INTENSITY, pa.float64(), nullable=True),
            pa.field(PRECURSOR_TOTAL_INTENSITY, pa.float64(), nullable=True),
            pa.field(ACQUISITION_DATETIME, pa.string(), nullable=True),
        ]
    )


def create_id_schema() -> pa.Schema:
    """Create and return the ID-only schema"""
    return pa.schema(
        [
            (SCAN, pa.string()),
            (MS_LEVEL, pa.int32()),
            (MZ_ARRAY, pa.list_(pa.float64())),
            (INTENSITY_ARRAY, pa.list_(pa.float64())),
        ]
    )


def batch_write_bruker_d(file_name: str, output_path: str, batch_size: int = 10000) -> str:
    """
    Batch processing and writing of Bruker .d files.

    Parameters
    ----------
    file_name : str
        Path to the Bruker .d directory
    output_path : str
        Path where the output parquet file will be saved
    batch_size : int
        Number of rows to process in each batch

    Returns
    -------
    str
        Path to the output parquet file
    """
    sql_filepath = f"{file_name}/analysis.tdf"
    schema = create_ms_schema()

    logger.info(f"Processing Bruker .d file: {file_name}")

    try:
        with sqlite3.connect(sql_filepath) as conn:
            # Retrieve acquisition datetime
            acquisition_date_time = conn.execute(
                "SELECT Value FROM GlobalMetadata WHERE key='AcquisitionDateTime'"
            ).fetchone()[0]

            # Check which optional columns exist
            columns = column_exists(conn, "frames")

            # Get allowed columns from the schema
            allowed_columns = {
                "Id": ("Id", SCAN),
                "MsMsType": (
                    "CASE WHEN MsMsType IN (8, 9) THEN 2 WHEN MsMsType = 0 THEN 1 ELSE NULL END",
                    MS_LEVEL,
                ),
                "NumPeaks": ("NumPeaks", NUM_PEAKS),
                "MaxIntensity": ("MaxIntensity", BASE_PEAK_INTENSITY),
                "SummedIntensities": ("SummedIntensities", SUMMED_PEAK_INTENSITY),
                "Time": ("Time", RETENTION_TIME),
                "Charge": ("Charge", CHARGE),
                "MonoisotopicMz": ("MonoisotopicMz", EXPERIMENTAL_MASS_TO_CHARGE),
            }

            # Construct safe column list
            safe_columns = []
            for schema_col_name, sql_expr in allowed_columns.items():
                if schema_col_name in columns or schema_col_name == "Id":
                    safe_columns.append(sql_expr[0])

            # Construct the query using safe columns
            query = f"""SELECT {', '.join(safe_columns)} FROM frames"""

            # Set up parquet writer with appropriate schema
            with pq.ParquetWriter(
                output_path, schema=schema, compression="gzip"
            ) as parquet_writer:
                # Stream data in batches
                for chunk in pd.read_sql_query(query, conn, chunksize=batch_size):
                    chunk[ACQUISITION_DATETIME] = acquisition_date_time
                    # Change column names to match the schema using allowed columns mapping
                    chunk.rename(
                        columns={v[0]: v[1] for v in allowed_columns.values()}, inplace=True
                    )
                    chunk[SCAN] = chunk[SCAN].astype(str)
                    for col in schema.names:
                        if col not in chunk.columns:
                            chunk[col] = None
                    batch_table = pa.Table.from_pandas(chunk, schema=schema)
                    parquet_writer.write_table(batch_table)

        logger.info(f"Successfully wrote Bruker data to {output_path}")
        return output_path

    except Exception as e:
        logger.error(f"Error processing Bruker .d file: {e}")
        raise


def column_exists(conn: sqlite3.Connection, table_name: str) -> Set[str]:
    """
    Fetch the existing columns in the specified SQLite table.

    Parameters
    ----------
    conn : sqlite3.Connection
        SQLite connection object
    table_name : str
        Name of the table to check

    Returns
    -------
    Set[str]
        Set of column names that exist in the table
    """
    table_info = pd.read_sql_query(f"PRAGMA table_info({table_name});", conn)
    return set(table_info["name"].tolist())


class BatchWritingConsumer:
    """
    A class to consume mass spectrometry data and write to a parquet file in batches.
    Processes mzML files using pyopenms streaming for memory efficiency.
    """

    def __init__(
        self,
        parquet_schema: pa.Schema,
        id_parquet_schema: pa.Schema,
        output_path: str,
        batch_size: int = 10000,
        ms2_file: bool = False,
        id_output_path: Optional[str] = None,
    ):
        """
        Initialize the BatchWritingConsumer.

        Parameters
        ----------
        parquet_schema : pa.Schema
            Schema for the main parquet output
        id_parquet_schema : pa.Schema
            Schema for the ID-only parquet output
        output_path : str
            Path for the main parquet output
        batch_size : int
            Number of rows to write in each batch
        ms2_file : bool
            Whether to generate ms2_file with the corresponding spectra
        id_output_path : Optional[str]
            Path for the ID-only parquet output
        """
        self.parquet_schema = parquet_schema
        self.id_parquet_schema = id_parquet_schema
        self.output_path = output_path
        self.id_output_path = id_output_path
        self.batch_size = batch_size
        self.ms2_file = ms2_file
        self.batch_data = []
        self.psm_parts = []
        self.parquet_writer = None
        self.id_parquet_writer = None
        self.acquisition_datetime = None

    def transform_mzml_spectrum(
        self,
        i: int,
        spectrum: oms.MSSpectrum,
        mzml_exp: oms.MSExperiment,
    ) -> None:
        """
        Process a spectrum and add it to the batch data.

        Parameters
        ----------
        i : int
            Index of the spectrum
        spectrum : oms.MSSpectrum
            Spectrum to process
        mzml_exp : oms.MSExperiment
            The experiment containing the spectrum
        """
        # Extract peaks data
        mz_array, intensity_array = spectrum.get_peaks()
        peak_count = len(mz_array)
        scan_id = extract_scan_id(spectrum)

        # Basic spectrum properties
        ms_level = spectrum.getMSLevel()
        rt = spectrum.getRT()

        # Skip unsupported MS levels
        if ms_level not in (1, 2):
            logger.debug(f"Skipping spectrum with unsupported MS level: {ms_level}")
            return

        # Calculate intensities if peaks exist
        if peak_count > 0:
            base_peak_intensity = float(np.max(intensity_array))
            total_intensity = float(np.sum(intensity_array))
        else:
            base_peak_intensity = None
            total_intensity = None

        # Build row data based on MS level
        if ms_level == 2 and spectrum.getPrecursors():
            # Process MS2 with precursors
            first_precursor_calculated = self._extract_first_precursor_data(spectrum, i, mzml_exp)

            # Extract spectrum ID for ID-only mode
            if self.ms2_file:
                self.psm_parts.append(
                    {
                        SCAN: scan_id,
                        MS_LEVEL: ms_level,
                        MZ_ARRAY: mz_array.tolist(),
                        INTENSITY_ARRAY: intensity_array.tolist(),
                    }
                )

            # Working only with the first precursor
            first_precursor = spectrum.getPrecursors()[0]
            charge_state = first_precursor.getCharge()
            exp_mz = first_precursor.getMZ()
            intensity = first_precursor.getIntensity()

            if intensity is None or intensity == 0.0 and first_precursor_calculated:
                # Use calculated precursor intensity if available
                intensity = first_precursor_calculated["intensity"]

            row_data = {
                SCAN: scan_id,
                MS_LEVEL: ms_level,
                NUM_PEAKS: peak_count,
                BASE_PEAK_INTENSITY: base_peak_intensity,
                SUMMED_PEAK_INTENSITY: total_intensity,
                RETENTION_TIME: float(rt),
                CHARGE: int(charge_state) if charge_state else None,
                EXPERIMENTAL_MASS_TO_CHARGE: float(exp_mz) if exp_mz else None,
                PRECURSOR_RT: (
                    first_precursor_calculated["rt"] if first_precursor_calculated else None
                ),
                INTENSITY: float(intensity) if intensity else None,
                PRECURSOR_TOTAL_INTENSITY: (
                    first_precursor_calculated["total_intensity"]
                    if first_precursor_calculated
                    else None
                ),
                ACQUISITION_DATETIME: (
                    str(self.acquisition_datetime) if self.acquisition_datetime else None
                ),
            }
        else:
            # Process MS1
            row_data = {
                SCAN: scan_id,
                MS_LEVEL: ms_level,
                CHARGE: None,
                NUM_PEAKS: peak_count,
                BASE_PEAK_INTENSITY: base_peak_intensity,
                SUMMED_PEAK_INTENSITY: total_intensity,
                RETENTION_TIME: float(rt),
                EXPERIMENTAL_MASS_TO_CHARGE: None,
                INTENSITY: None,
                PRECURSOR_RT: None,
                PRECURSOR_TOTAL_INTENSITY: None,
                ACQUISITION_DATETIME: (
                    str(self.acquisition_datetime) if self.acquisition_datetime else None
                ),
            }

        self.batch_data.append(row_data)

        # Write batch when it reaches specified size
        if len(self.batch_data) >= self.batch_size:
            self._write_batch()

    def _extract_first_precursor_data(
        self, spectrum: oms.MSSpectrum, i: int, mzml_exp: oms.MSExperiment
    ) -> Dict:
        """
        Extract precursor information from the first precursor when it is not annotated in the MS2.

        Parameters
        ----------
        spectrum : oms.MSSpectrum
            Spectrum to extract precursors from
        i : int
            Index of the spectrum
        mzml_exp : oms.MSExperiment
            The experiment containing the spectrum

        Returns
        -------
        Dict
            First precursor data
        """
        first_precursor = {}

        if spectrum.getPrecursors():
            index = 0
            precursor = spectrum.getPrecursors()[index]

            # Get precursor spectrum details if available
            precursor_spectrum_index = mzml_exp.getPrecursorSpectrum(i)
            precursor_rt = None
            precursor_intensity = None
            total_intensity = None

            if precursor_spectrum_index >= 0:
                precursor_spectrum = mzml_exp.getSpectrum(precursor_spectrum_index)
                precursor_rt = precursor_spectrum.getRT()

                # Calculate purity metrics
                try:
                    purity = oms.PrecursorPurity().computePrecursorPurity(
                        precursor_spectrum, precursor, 100, True
                    )
                    precursor_intensity = purity.target_intensity
                    total_intensity = purity.total_intensity
                except Exception as e:
                    logger.debug(f"Could not compute precursor purity: {e}")

            first_precursor = {
                "index": index,
                "rt": precursor_rt,
                "intensity": precursor_intensity,
                "total_intensity": total_intensity,
            }
        else:
            logger.debug(f"No precursors found for spectrum: {spectrum.getNativeID()}")

        return first_precursor

    def _write_batch(self) -> None:
        """
        Write accumulated batch data efficiently using PyArrow's streaming writer.
        """
        try:
            # Write main data batch if exists
            if self.batch_data:
                # Initialize writer lazily if not already created
                if self.parquet_writer is None:
                    self.parquet_writer = pq.ParquetWriter(
                        where=self.output_path, schema=self.parquet_schema, compression="gzip"
                    )

                # Create a RecordBatch directly from the current batch
                batch = pa.RecordBatch.from_pylist(self.batch_data, schema=self.parquet_schema)

                # Write the batch directly
                self.parquet_writer.write_batch(batch)
                self.batch_data = []

            # Write ID-only data if enabled and data exists
            if self.ms2_file and self.psm_parts:
                # Initialize writer lazily
                if self.id_parquet_writer is None:
                    self.id_parquet_writer = pq.ParquetWriter(
                        where=self.id_output_path,
                        schema=self.id_parquet_schema,
                        compression="gzip",
                    )

                # Create and write batch
                batch = pa.RecordBatch.from_pylist(self.psm_parts, schema=self.id_parquet_schema)
                self.id_parquet_writer.write_batch(batch)
                self.psm_parts = []

        except Exception as e:
            logger.error(f"Error during batch writing: {e}, file path: {self.output_path}")
            raise

    def finalize(self) -> None:
        """
        Finalize the writing process by writing any remaining data and closing writers.
        """
        # Write any remaining data
        if self.batch_data or self.psm_parts:
            self._write_batch()

        # Close writers if they exist
        if self.parquet_writer:
            self.parquet_writer.close()
            self.parquet_writer = None

        if self.id_parquet_writer:
            self.id_parquet_writer.close()
            self.id_parquet_writer = None


def batch_write_mzml_streaming(
    file_name: str,
    output_path: str,
    ms2_file: bool = False,
    id_output_path: Optional[str] = None,
    batch_size: int = 10000,
) -> Optional[str]:
    """
    Process an mzML file using streaming and write to parquet.

    Parameters
    ----------
    file_name : str
        Path to the mzML file
    output_path : str
        Path for the main output
    ms2_file : bool
        Whether to generate ID-only output
    id_output_path : Optional[str]
        Path for the ID-only output
    batch_size : int
        Number of rows to process in each batch

    Returns
    -------
    Optional[str]
        Path to the output file or None if processing failed
    """
    logger.info(f"Processing mzML file: {file_name}")

    try:
        # Create schemas
        parquet_schema = create_ms_schema()
        id_parquet_schema = create_id_schema()

        # Set up mzML processing
        mzml_handler = oms.MzMLFile()
        mzml_exp = oms.MSExperiment()
        mzml_handler.load(file_name, mzml_exp)

        # Create batch writer
        batch_writer = BatchWritingConsumer(
            parquet_schema=parquet_schema,
            id_parquet_schema=id_parquet_schema,
            output_path=output_path,
            ms2_file=ms2_file,
            id_output_path=id_output_path,
            batch_size=batch_size,
        )

        # Get acquisition datetime if available
        if mzml_exp.getMetaValue("acquisition_date_time"):
            batch_writer.acquisition_datetime = mzml_exp.getMetaValue("acquisition_date_time")
        else:
            acquisition = mzml_exp.getDateTime().get()
            if acquisition:
                batch_writer.acquisition_datetime = acquisition

        # Process each spectrum
        for i, spec in enumerate(mzml_exp):
            batch_writer.transform_mzml_spectrum(i, spec, mzml_exp)

        # Finalize processing
        batch_writer.finalize()
        logger.info(f"Successfully processed mzML file to {output_path}")
        return output_path

    except Exception as e:
        logger.error(f"Error during mzML streaming: {e}, file path: {file_name}")
        return None


def resolve_ms_path(ms_path: str) -> str:
    """
    Resolve mass spectrometry file path with improved candidate search.

    Parameters
    ----------
    ms_path : str
        Path to resolve

    Returns
    -------
    str
        Resolved file path

    Raises
    ------
    FileNotFoundError
        If no matching file is found
    """
    path_obj = Path(ms_path)
    if path_obj.exists():
        return str(path_obj)

    # Look for files with matching stem and valid extensions
    candidates = list(path_obj.parent.glob(f"{path_obj.stem}*"))
    valid_extensions = {".d", ".mzml", ".mzML"}
    valid_candidates = [
        str(c.resolve()) for c in candidates if c.suffix.lower() in valid_extensions
    ]

    if len(valid_candidates) == 1:
        logger.info(f"Resolved {ms_path} to {valid_candidates[0]}")
        return valid_candidates[0]
    elif len(valid_candidates) > 1:
        logger.warning(f"Multiple candidates found for {ms_path}: {valid_candidates}")

    raise FileNotFoundError(f"No unique file found for {ms_path}")


@click.command("mzmlstats")
@click.option(
    "--ms_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to mass spectrometry file (.mzML or .d)",
)
@click.option(
    "--ms2_file", is_flag=True, help="Generate a parquet with the spectrum id and the peaks"
)
@click.option(
    "--feature_detection",
    is_flag=True,
    help="Run feature detection on MS1 and get the features file",
)
@click.option(
    "--batch_size", type=int, default=10000, help="Number of rows to write in each batch"
)
@click.pass_context
def mzml_statistics(
    ctx,
    ms_path: str,
    ms2_file: bool = False,
    feature_detection: bool = False,
    batch_size: int = 10000,
) -> None:
    """
    Parse mass spectrometry data files (.mzML or Bruker .d formats) to extract
    and compile statistics about the spectra.

    Example usage:
    quantmsutilsc mzmlstats --ms_path "path/to/file.mzML" --ms2_file --batch_size 5000
    """
    try:
        # Resolve the file path
        ms_path = resolve_ms_path(ms_path)
        path_obj = Path(ms_path)

        # Prepare output paths
        output_path = str(path_obj.with_name(f"{path_obj.stem}_ms_info.parquet"))
        id_output_path = str(path_obj.with_name(f"{path_obj.stem}_ms2_info.parquet"))
        feature_output_path = str(path_obj.with_name(f"{path_obj.stem}_ms1_feature_info.parquet"))

        # Process based on file type
        if path_obj.suffix.lower() == ".d":
            batch_write_bruker_d(file_name=ms_path, output_path=output_path, batch_size=batch_size)
        if path_obj.suffix.lower() in [".mzml"]:
            batch_write_mzml_streaming(
                file_name=ms_path,
                output_path=output_path,
                ms2_file=ms2_file,
                id_output_path=id_output_path,
                batch_size=batch_size,
            )
            if feature_detection:
                feature_detector = MS1FeatureDetector()
                feature_detector.process_file(input_file=ms_path, output_file=feature_output_path)
            logger.info("The file {} has been processed".format(ms_path))
        else:
            raise ValueError(f"Unsupported file type: {path_obj.suffix}")

        logger.info(f"Successfully processed mass spectrometry file: {ms_path}")

    except Exception as e:
        logger.error(f"Error processing file {ms_path}: {e}")
        raise click.Abort()
