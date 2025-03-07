import logging
import re
import sqlite3
from pathlib import Path
from typing import Optional, Set, Union

import click
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyopenms as oms
from pyarrow import Schema

from quantmsutils.utils.constants import (
    CHARGE,
    SCAN,
    MS_LEVEL,
    NUM_PEAKS,
    BASE_PEAK_INTENSITY,
    SUMMED_PEAK_INTENSITY,
    RETENTION_TIME,
    EXPERIMENTAL_MASS_TO_CHARGE,
    ACQUISITION_DATETIME,
    MZ_ARRAY,
    INTENSITY_ARRAY,
    MONOISOTOPIC_MZ,
    MAX_INTENSITY,
    PRECURSORS,
    INTENSITY,
)

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def batch_write_bruker_d(
    file_name: str, output_path: str, batch_size: int = 10000, schema: Schema = None
) -> str:
    """
    Batch processing and writing of Bruker .d files.
    """
    sql_filepath = f"{file_name}/analysis.tdf"

    with sqlite3.connect(sql_filepath) as conn:
        # Retrieve acquisition datetime
        acquisition_date_time = conn.execute(
            "SELECT Value FROM GlobalMetadata WHERE key='AcquisitionDateTime'"
        ).fetchone()[0]

        # Check which optional columns exist
        columns = column_exists(conn, "frames")

        # Get allowed columns from the schema
        allowed_columns = {
            "Id": "Id",
            "MsMsType": "CASE WHEN MsMsType IN (8, 9) THEN 2 WHEN MsMsType = 0 THEN 1 ELSE NULL END",
            "NumPeaks": "NumPeaks",
            "MaxIntensity": "MaxIntensity",
            "SummedIntensities": "SummedIntensities",
            "Time": "Time",
            "Charge": "Charge",
            "MonoisotopicMz": "MonoisotopicMz",
        }

        # Construct safe column list
        safe_columns = []
        column_mapping = {}
        for schema_col_name, sql_expr in allowed_columns.items():
            if schema_col_name in columns or schema_col_name == "Id":
                safe_columns.append(sql_expr)
                column_mapping[schema_col_name] = sql_expr

        # Construct the query using parameterized safe columns
        query = f"""SELECT {', '.join(safe_columns)} FROM frames"""

        schema = pa.schema(
            [
                pa.field(SCAN, pa.int32(), nullable=False),
                pa.field(MS_LEVEL, pa.int32(), nullable=True),
                pa.field(NUM_PEAKS, pa.int32(), nullable=True),
                pa.field(MAX_INTENSITY, pa.float64(), nullable=True),
                pa.field(SUMMED_PEAK_INTENSITY, pa.float64(), nullable=True),
                pa.field(RETENTION_TIME, pa.float64(), nullable=True),
                pa.field(CHARGE, pa.int32(), nullable=True),
                pa.field(MONOISOTOPIC_MZ, pa.float64(), nullable=True),
                pa.field(INTENSITY, pa.float64(), nullable=True),
                pa.field(
                    PRECURSORS,
                    pa.list_(
                        pa.struct(
                            [
                                ("charge", pa.int32()),
                                ("mz", pa.float64()),
                                ("intensity", pa.float64()),
                                ("rt", pa.float64()),
                                ("index", pa.int32()),
                            ]
                        )
                    ),
                    nullable=True,
                ),
                pa.field(ACQUISITION_DATETIME, pa.string(), nullable=True),
            ]
        )

        # Set up parquet writer
        parquet_writer = pq.ParquetWriter(output_path, schema=schema, compression="gzip")

        try:
            # Stream data in batches
            for chunk in pd.read_sql_query(query, conn, chunksize=batch_size):
                chunk["AcquisitionDateTime"] = acquisition_date_time
                for col in schema.names:
                    if col not in chunk.columns:
                        chunk[col] = None
                batch_table = pa.Table.from_pandas(chunk, schema=schema)
                parquet_writer.write_table(batch_table)

        finally:
            parquet_writer.close()

    return output_path


def column_exists(conn, table_name: str) -> Set[str]:
    """
    Fetch the existing columns in the specified SQLite table.
    """
    table_info = pd.read_sql_query(f"PRAGMA table_info({table_name});", conn)
    return set(table_info["name"].tolist())


class BatchWritingConsumer:
    """
    A class to consume mass spectrometry data and write to a parquet file in batches from mzML files using
    pyopenms streaming.
    """

    def __init__(
        self,
        parquet_schema: pa.Schema,
        id_parquet_schema: pa.Schema,
        output_path,
        batch_size=10000,
        id_only=False,
        id_output_path=None,
    ):
        self.parquet_schema = parquet_schema
        self.id_parquet_schema = id_parquet_schema
        self.output_path = output_path
        self.id_output_path = id_output_path
        self.batch_size = batch_size
        self.id_only = id_only
        self.batch_data = []
        self.psm_parts = []
        self.parquet_writer = None
        self.id_parquet_writer = None
        self.acquisition_datetime = None
        self.scan_pattern = re.compile(r"[scan|spectrum]=(\d+)")

    def transform_mzml_spectrum(
        self,
        i: int,
        spectrum: oms.MSSpectrum,
        mzml_exp: oms.MSExperiment,
        mzml_handler: oms.MzMLFile,
    ):
        """
        Process a given spectrum to extract relevant data and append it to the batch data.

        This method extracts peak information, MS level, retention time, and other
        relevant attributes from the spectrum. It organizes the data into a dictionary
        and appends it to the batch data. If the MS level is 2 and `id_only` is True,
        it also extracts precursor information and appends it to `psm_parts`.
        When the batch data reaches the specified batch size, it triggers a batch write.

        Parameters
        ----------
        spectrum : pyopenms.MSSpectrum
            The spectrum object to process.
        """

        mz_array, intensity_array = (
            spectrum.get_peaks()
        )  # Split the peaks into mz and intensity arrays
        peak_per_ms = len(mz_array)  # Get the number of peaks in the spectrum
        base_peak_intensity = (
            float(np.max(intensity_array)) if peak_per_ms > 0 else None
        )  # Get the maximum intensity
        total_intensity = (
            float(np.sum(intensity_array)) if peak_per_ms > 0 else None
        )  # Get the total intensity # TODO: Review by @timo and @julianus
        ms_level = spectrum.getMSLevel()  # Get the MS level of the spectrum
        rt = spectrum.getRT()  # Get the retention time of the spectrum

        if ms_level == 2 and spectrum.getPrecursors():
            precursor = spectrum.getPrecursors()[0]  # Get the first precursor
            charge_state = precursor.getCharge()  # Charge of first precursor
            exp_mz = precursor.getMZ()  # Experimental mass to charge ratio of first precursor
            intensity = precursor.getIntensity()  # Intensity of first precursor
            precursors = []
            # capture if more than one precursor
            index = 0
            for precursor in spectrum.getPrecursors():

                logging.info(
                    "Precursor charge: %s, precursor mz: %s",
                    precursor.getCharge(),
                    precursor.getMZ(),
                )
                charge_state = precursor.getCharge()  # TODO: Review by @timo and @julianus
                exp_mz = precursor.getMZ()  # TODO: Review by @timo and @julianus
                intensity = precursor.getIntensity()  # TODO: Review by @timo and @julianus
                precursor_spectrum_index = mzml_exp.getPrecursorSpectrum(i)
                precursor_rt = None
                total_intensity = None
                precursor_intensity = None
                precursor_intensity_in_isolation_window = None
                if precursor_spectrum_index >= 0:
                    precursor_spectrum = mzml_exp.getSpectrum(precursor_spectrum_index)
                    precursor_rt = precursor_spectrum.getRT()
                    purity = oms.PrecursorPurity().computePrecursorPurity(
                        precursor_spectrum, precursor, 20.0, True
                    )
                    precursor_intensity = purity.target_intensity
                    precursor_intensity_in_isolation_window = purity.total_intensity
                precursors.append(
                    {
                        "index": index,
                        "charge": charge_state,
                        "mz": exp_mz,
                        "rt": precursor_rt,
                        "intensity": precursor_intensity,
                        "total_intensity": total_intensity,
                        "intensity_isolation_window": precursor_intensity_in_isolation_window,
                    }
                )
                index += 1

            if self.id_only:
                scan_id = self.scan_pattern.findall(spectrum.getNativeID())[0]
                self.psm_parts.append(
                    {
                        SCAN: scan_id,
                        MS_LEVEL: int(ms_level),
                        MZ_ARRAY: mz_array.tolist(),
                        INTENSITY_ARRAY: intensity_array.tolist(),
                    }
                )

            row_data = {
                SCAN: spectrum.getNativeID(),
                MS_LEVEL: int(ms_level),
                CHARGE: int(charge_state) if charge_state is not None else None,
                NUM_PEAKS: int(peak_per_ms),
                BASE_PEAK_INTENSITY: (
                    float(base_peak_intensity) if base_peak_intensity is not None else None
                ),
                SUMMED_PEAK_INTENSITY: (
                    float(total_intensity) if total_intensity is not None else None
                ),
                RETENTION_TIME: float(rt),
                EXPERIMENTAL_MASS_TO_CHARGE: float(exp_mz) if exp_mz is not None else None,
                INTENSITY: float(intensity) if intensity is not None else None,
                ACQUISITION_DATETIME: str(self.acquisition_datetime),
                PRECURSORS: precursors,
            }
        elif ms_level == 1:
            row_data = {
                SCAN: spectrum.getNativeID(),
                MS_LEVEL: int(ms_level),
                CHARGE: None,
                NUM_PEAKS: int(peak_per_ms),
                BASE_PEAK_INTENSITY: (
                    float(base_peak_intensity) if base_peak_intensity is not None else None
                ),
                SUMMED_PEAK_INTENSITY: (
                    float(total_intensity) if total_intensity is not None else None
                ),
                RETENTION_TIME: float(rt),
                EXPERIMENTAL_MASS_TO_CHARGE: None,
                ACQUISITION_DATETIME: str(self.acquisition_datetime),
                PRECURSORS: None,
            }
        else:
            logger.info(
                "Skipping spectrum with MS level %s, MS not supported, or precursors not ",
                ms_level,
            )
            return

        self.batch_data.append(row_data)

        # Write batch when it reaches specified size
        if len(self.batch_data) >= self.batch_size:
            self._write_batch()

    def _write_batch(self):
        """
        Write accumulated batch data more efficiently using PyArrow's streaming writer.

        Improvements:
        - Directly stream data without creating a full in-memory table
        - Reduce memory overhead for large datasets
        - More efficient batch processing
        """
        try:

            if self.batch_data:

                # Initialize writers lazily if not already created
                if self.parquet_writer is None:
                    self.parquet_writer = pq.ParquetWriter(
                        where=self.output_path, schema=self.parquet_schema, compression="gzip"
                    )

                # Create a Table directly from the current batch
                batch = pa.RecordBatch.from_pylist(self.batch_data, schema=self.parquet_schema)

                # Write the batch directly
                self.parquet_writer.write_batch(batch)

                # Clear the batch data
                self.batch_data = []

            if self.id_only and self.psm_parts:

                # Initialize writers lazily if not already created
                if self.id_parquet_writer is None:
                    self.id_parquet_writer = pq.ParquetWriter(
                        where=self.id_output_path,
                        schema=self.id_parquet_schema,
                        compression="gzip",
                    )

                # Create a Table directly from the current batch
                batch = pa.RecordBatch.from_pylist(self.psm_parts, schema=self.id_parquet_schema)
                # Write the batch directly
                self.id_parquet_writer.write_batch(batch)
                self.psm_parts = []

        except Exception as e:
            logger.error(f"Error during batch writing: {e}, file path: {self.output_path}")
            raise

    def finalize(self):
        """
        Finalize the writing process.
        :return:
        """
        if self.batch_data or self.psm_parts:
            self._write_batch()

        if self.parquet_writer:
            self.parquet_writer.close()

        if self.id_parquet_writer:
            self.id_parquet_writer.close()


def batch_write_mzml_streaming(
    file_name: str,
    parquet_schema: Schema,
    id_parquet_schema: Schema,
    output_path: str,
    id_only: bool,
    id_output_path: str,
    batch_size: int = 10000,
) -> Union[str, None]:

    try:
        mzml_handler = oms.MzMLFile()
        mzml_exp = oms.MSExperiment()
        mzml_handler.load(file_name, mzml_exp)
        batch_writer = BatchWritingConsumer(
            parquet_schema=parquet_schema,
            id_parquet_schema=id_parquet_schema,
            output_path=output_path,
            id_only=id_only,
            id_output_path=id_output_path,
            batch_size=batch_size,
        )

        for i, spec in enumerate(mzml_exp):
            batch_writer.transform_mzml_spectrum(i, spec, mzml_exp, mzml_handler)

        return output_path
    except Exception as e:
        logger.error(f"Error during streaming: {e}, file path: {file_name}")
        return None


@click.command("mzmlstats")
@click.option("--ms_path", type=click.Path(exists=True), required=True)
@click.option("--id_only", is_flag=True, help="Generate a csv with the spectrum id and the peaks")
@click.option(
    "--batch_size", type=int, default=10000, help="Number of rows to write in each batch"
)
@click.pass_context
def mzml_statistics(ctx, ms_path: str, id_only: bool = False, batch_size: int = 10000) -> None:
    """
    The mzml_statistics function parses mass spectrometry data files, either in
    .mzML or Bruker .d formats, to extract and compile a set of statistics about the spectra contained within.
    It supports generating detailed or ID-only CSV files based on the spectra data.

    # Command line usage example
    quantmsutilsc mzmlstats --ms_path "path/to/file.mzML"

    :param ctx: Click context

    :param ms_path: A string specifying the path to the mass spectrometry file.
    :param id_only: A boolean flag that, when set to True, generates a CSV file containing only the spectrum ID and
    peaks data for MS level 2 spectra.
    :param batch_size: An integer specifying the number of rows to write in each batch.

    """
    schema = pa.schema(
        [
            pa.field(SCAN, pa.string(), nullable=True),
            pa.field(MS_LEVEL, pa.int32(), nullable=True),
            pa.field(CHARGE, pa.int32(), nullable=True),
            pa.field(NUM_PEAKS, pa.int32(), nullable=True),
            pa.field(BASE_PEAK_INTENSITY, pa.float64(), nullable=True),
            pa.field(SUMMED_PEAK_INTENSITY, pa.float64(), nullable=True),
            pa.field(RETENTION_TIME, pa.float64(), nullable=True),
            pa.field(EXPERIMENTAL_MASS_TO_CHARGE, pa.float64(), nullable=True),
            pa.field(INTENSITY, pa.float64(), nullable=True),
            pa.field(
                PRECURSORS,
                pa.list_(
                    pa.struct(
                        [
                            ("charge", pa.int32()),
                            ("mz", pa.float64()),
                            ("intensity", pa.float64()),
                            ("rt", pa.float64()),
                            ("total_intensity", pa.float64()),
                            ("intensity_isolation_window", pa.float64()),
                            ("index", pa.int32()),
                        ]
                    )
                ),
                nullable=True,
            ),
            pa.field(ACQUISITION_DATETIME, pa.string(), nullable=True),
        ]
    )

    id_schema = pa.schema(
        [
            (SCAN, pa.string()),
            (MS_LEVEL, pa.int32()),
            (MZ_ARRAY, pa.list_(pa.float64())),
            (INTENSITY_ARRAY, pa.list_(pa.float64())),
        ]
    )

    ms_path = _resolve_ms_path(ms_path)
    output_path = str(Path(ms_path).with_name(f"{Path(ms_path).stem}_ms_info.parquet"))
    id_output_path = str(Path(ms_path).with_name(f"{Path(ms_path).stem}_spectrum_df.parquet"))

    if Path(ms_path).suffix == ".d":
        batch_write_bruker_d(file_name=ms_path, output_path=output_path, batch_size=batch_size)
    elif Path(ms_path).suffix.lower() in [".mzml"]:
        batch_write_mzml_streaming(
            file_name=ms_path,
            parquet_schema=schema,
            id_parquet_schema=id_schema,
            output_path=output_path,
            id_only=id_only,
            id_output_path=id_output_path,
            batch_size=batch_size,
        )
    else:
        raise RuntimeError(f"Unsupported file type: {ms_path}")

    logging.info("Successfully processed mass spectrometry file: %s", ms_path)


def _resolve_ms_path(ms_path: str) -> str:
    """
    Resolve mass spectrometry file path with improved candidate search.
    """
    path_obj = Path(ms_path)
    if path_obj.exists():
        return str(path_obj)

    candidates = list(path_obj.parent.glob(f"{path_obj.stem}*"))
    valid_extensions = {".d", ".mzml", ".mzML"}
    candidates = [str(c.resolve()) for c in candidates if c.suffix.lower() in valid_extensions]

    if len(candidates) == 1:
        return candidates[0]

    raise FileNotFoundError(f"No unique file found for {ms_path}")
