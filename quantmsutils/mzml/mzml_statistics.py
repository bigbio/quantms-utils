import re
import sqlite3
from pathlib import Path
from typing import Optional, List

import click
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyopenms import MzMLFile


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
    ):
        self.parquet_schema = parquet_schema
        self.id_parquet_schema = id_parquet_schema
        self.output_path = output_path
        self.batch_size = batch_size
        self.id_only = id_only
        self.batch_data = []
        self.psm_parts = []
        self.parquet_writer = None
        self.id_parquet_writer = None
        self.acquisition_datetime = None
        self.scan_pattern = re.compile(r"[scan|spectrum]=(\d+)")

    def setExperimentalSettings(self, settings):
        self.acquisition_datetime = settings.getDateTime().get()

    def setExpectedSize(self, a, b):
        pass

    def consumeChromatogram(self, chromatogram):
        pass

    def consumeSpectrum(self, spectrum):
        """
        Consume spectrum data and write to parquet file.
        :param spectrum: spectrum data.
        :return: None
        """

        peaks = spectrum.get_peaks()
        mz_array, intensity_array = peaks[0], peaks[1]
        peak_per_ms = len(mz_array)
        base_peak_intensity = float(np.max(intensity_array)) if peak_per_ms > 0 else None
        total_intensity = float(np.sum(intensity_array)) if peak_per_ms > 0 else None
        ms_level = spectrum.getMSLevel()
        rt = spectrum.getRT()

        if ms_level == 2:
            precursor = spectrum.getPrecursors()[0]
            charge_state = precursor.getCharge()
            exp_mz = precursor.getMZ()

            if self.id_only:
                scan_id = self.scan_pattern.findall(spectrum.getNativeID())[0]
                self.psm_parts.append(
                    [
                        {
                            "scan": scan_id,
                            "ms_level": ms_level,
                            "mz": mz_array,
                            "intensity": intensity_array,
                        }
                    ]
                )

            row_data = {
                "SpectrumID": spectrum.getNativeID(),
                "MSLevel": float(ms_level),
                "Charge": float(charge_state) if charge_state is not None else None,
                "MS_peaks": float(peak_per_ms),
                "Base_Peak_Intensity": (
                    float(base_peak_intensity) if base_peak_intensity is not None else None
                ),
                "Summed_Peak_Intensities": (
                    float(total_intensity) if total_intensity is not None else None
                ),
                "Retention_Time": float(rt),
                "Exp_Mass_To_Charge": float(exp_mz) if exp_mz is not None else None,
                "AcquisitionDateTime": str(self.acquisition_datetime),
            }
        elif ms_level == 1:
            row_data = {
                "SpectrumID": spectrum.getNativeID(),
                "MSLevel": float(ms_level),
                "Charge": None,
                "MS_peaks": float(peak_per_ms),
                "Base_Peak_Intensity": (
                    float(base_peak_intensity) if base_peak_intensity is not None else None
                ),
                "Summed_Peak_Intensities": (
                    float(total_intensity) if total_intensity is not None else None
                ),
                "Retention_Time": float(rt),
                "Exp_Mass_To_Charge": None,
                "AcquisitionDateTime": str(self.acquisition_datetime),
            }
        else:
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
            # If no data, return early
            if not self.batch_data:
                return

            # Initialize writers lazily if not already created
            if self.parquet_writer is None:
                self.parquet_writer = pq.ParquetWriter(
                    where=self.output_path, schema=self.parquet_schema, compression="gzip"
                )

            # Create a RecordBatch directly from the current batch
            batch = pa.RecordBatch.from_pylist(self.batch_data, schema=self.parquet_schema)

            # Write the batch directly
            self.parquet_writer.write_batch(batch)

            # Clear the batch data
            self.batch_data = []

            # Handle ID-only data if applicable
            if self.id_only and self.psm_parts:
                # Similar approach for spectrum ID data
                if self.id_parquet_writer is None:
                    self.id_parquet_writer = pq.ParquetWriter(
                        where=f"{Path(self.output_path).stem}_spectrum_df.parquet",
                        schema=self.id_parquet_schema,
                        compression="gzip",
                    )

                id_batch = pa.RecordBatch.from_pylist(
                    self.psm_parts, schema=self.id_parquet_schema
                )
                self.id_parquet_writer.write_batch(id_batch)
                self.psm_parts = []

        except Exception as e:
            print(f"Error during batch writing: {e}")
            raise

    def finalize(self):
        """
        Finalize the writing process.
        :return:
        """
        if self.batch_data:
            self._write_batch()

        # Write spectrum data if id_only
        if self.id_only and self.psm_parts:
            self._write_batch()

        if self.parquet_writer:
            self.parquet_writer.close()

        if self.id_parquet_writer:
            self.id_parquet_writer.close()


def column_exists(conn, table_name: str) -> List[str]:
    """
    Fetch the existing columns in the specified SQLite table.
    """
    table_info = pd.read_sql_query(f"PRAGMA table_info({table_name});", conn)
    return set(table_info["name"].tolist())


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
            pa.field("SpectrumID", pa.string(), nullable=True),
            pa.field("MSLevel", pa.float64(), nullable=True),
            pa.field("Charge", pa.float64(), nullable=True),
            pa.field("MS_peaks", pa.float64(), nullable=True),
            pa.field("Base_Peak_Intensity", pa.float64(), nullable=True),
            pa.field("Summed_Peak_Intensities", pa.float64(), nullable=True),
            pa.field("Retention_Time", pa.float64(), nullable=True),
            pa.field("Exp_Mass_To_Charge", pa.float64(), nullable=True),
            pa.field("AcquisitionDateTime", pa.string(), nullable=True),
        ]
    )

    id_schema = pa.schema(
        [
            ("scan", pa.string()),
            ("ms_level", pa.int32()),
            ("mz", pa.list_(pa.float64())),
            ("intensity", pa.list_(pa.float64())),
        ]
    )

    def batch_write_mzml_streaming(
        file_name: str,
        parquet_schema: pa.Schema,
        output_path: str,
        id_parquet_schema: pa.Schema,
        id_only: bool = False,
        batch_size: int = 10000,
    ) -> Optional[str]:
        """
        Parse mzML file in a streaming manner and write to Parquet.
        """
        consumer = BatchWritingConsumer(
            parquet_schema=parquet_schema,
            output_path=output_path,
            batch_size=batch_size,
            id_only=id_only,
            id_parquet_schema=id_parquet_schema,
        )
        try:
            MzMLFile().transform(file_name.encode(), consumer)
            consumer.finalize()
            return output_path
        except Exception as e:
            print(f"Error during streaming: {e}")
            return None

    def batch_write_bruker_d(file_name: str, output_path: str, batch_size: int = 10000) -> str:
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
                    pa.field("Id", pa.int32(), nullable=False),
                    pa.field("MsMsType", pa.int32(), nullable=True),
                    pa.field("NumPeaks", pa.int32(), nullable=True),
                    pa.field("MaxIntensity", pa.float64(), nullable=True),
                    pa.field("SummedIntensities", pa.float64(), nullable=True),
                    pa.field("Time", pa.float64(), nullable=True),
                    pa.field("Charge", pa.int32(), nullable=True),
                    pa.field("MonoisotopicMz", pa.float64(), nullable=True),
                    pa.field("AcquisitionDateTime", pa.string(), nullable=True),
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

    # Resolve file path
    ms_path = _resolve_ms_path(ms_path)
    output_path = f"{Path(ms_path).stem}_ms_info.parquet"

    # Choose processing method based on file type
    if Path(ms_path).suffix == ".d":
        batch_write_bruker_d(file_name=ms_path, output_path=output_path, batch_size=batch_size)
    elif Path(ms_path).suffix.lower() in [".mzml"]:
        batch_write_mzml_streaming(
            file_name=ms_path,
            parquet_schema=schema,
            id_parquet_schema=id_schema,
            output_path=output_path,
            id_only=id_only,
            batch_size=batch_size,
        )
    else:
        raise RuntimeError(f"Unsupported file type: {ms_path}")


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
