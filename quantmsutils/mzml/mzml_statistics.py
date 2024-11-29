import re
import sqlite3
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyopenms import MSExperiment, MzMLFile


@click.command("mzmlstats")
@click.option("--ms_path", type=click.Path(exists=True), required=True)
@click.option(
    "--id_only", is_flag=True,
    help="Generate a csv with the spectrum id and the peaks"
)
@click.option(
    "--batch_size",
    type=int,
    default=10000,
    help="Number of rows to write in each batch"
)
@click.pass_context
def mzml_statistics(
    ctx,
    ms_path: str,
    id_only: bool = False,
    batch_size: int = 10000
) -> None:
    """
    The mzml_statistics function parses mass spectrometry data files, either in
    .mzML or Bruker .d formats, to extract and compile a set of statistics about the spectra contained within.
    It supports generating detailed or ID-only CSV files based on the spectra data.

    # Command line usage example
    python script_name.py mzml_statistics --ms_path "path/to/file.mzML"

    :param ctx: Click context
    :param ms_path: A string specifying the path to the mass spectrometry file.
    :param id_only: A boolean flag that, when set to True, generates a CSV file containing only the spectrum ID and
    peaks data for MS level 2 spectra.

    """
    schema = pa.schema([
        pa.field("SpectrumID", pa.string(), nullable=True),
        pa.field("MSLevel", pa.float64(), nullable=True),
        pa.field("Charge", pa.float64(), nullable=True),
        pa.field("MS_peaks", pa.float64(), nullable=True),
        pa.field("Base_Peak_Intensity", pa.float64(), nullable=True),
        pa.field("Summed_Peak_Intensities", pa.float64(), nullable=True),
        pa.field("Retention_Time", pa.float64(), nullable=True),
        pa.field("Exp_Mass_To_Charge", pa.float64(), nullable=True),
        pa.field("AcquisitionDateTime", pa.string(), nullable=True),
    ])

    def batch_write_mzml(
        file_name: str,
        parquet_schema: pa.Schema,
        output_path: str,
        id_only: bool = False,
        batch_size: int = 10000
    ) -> Optional[str]:
        """
        Parse mzML file and return a pandas DataFrame with the information. If id_only is True, it will also save a csv.
        @param file_name: The file name of the mzML file
        @param id_only: If True, it will save a csv with the spectrum id, mz and intensity
        @return: A pandas DataFrame with the information of the mzML file
        """
        exp = MSExperiment()
        MzMLFile().load(file_name, exp)
        acquisition_datetime = exp.getDateTime().get()
        scan_pattern = re.compile(r"[scan|spectrum]=(\d+)")

        # Prepare parquet writer
        parquet_writer = None
        batch_data = []
        psm_parts = []

        try:
            for spectrum in exp:
                # Peak extraction
                peaks = spectrum.get_peaks()
                mz_array, intensity_array = peaks[0], peaks[1]

                # Compute peaks efficiently
                peak_per_ms = len(mz_array)
                base_peak_intensity = float(np.max(intensity_array)) if peak_per_ms > 0 else None
                total_intensity = float(np.sum(intensity_array)) if peak_per_ms > 0 else None

                # Metadata extraction
                ms_level = spectrum.getMSLevel()
                rt = spectrum.getRT()

                if ms_level == 2:
                    precursor = spectrum.getPrecursors()[0]
                    charge_state = precursor.getCharge()
                    exp_mz = precursor.getMZ()

                    if id_only:
                        scan_id = scan_pattern.findall(spectrum.getNativeID())[0]
                        psm_parts.append([
                            scan_id, ms_level,
                            mz_array.tolist(),
                            intensity_array.tolist()
                        ])

                    row_data = {
                        "SpectrumID": spectrum.getNativeID(),
                        "MSLevel": float(ms_level),
                        "Charge": float(charge_state) if charge_state is not None else None,
                        "MS_peaks": float(peak_per_ms),
                        "Base_Peak_Intensity": float(base_peak_intensity) if base_peak_intensity is not None else None,
                        "Summed_Peak_Intensities": float(total_intensity) if total_intensity is not None else None,
                        "Retention_Time": float(rt),
                        "Exp_Mass_To_Charge": float(exp_mz) if exp_mz is not None else None,
                        "AcquisitionDateTime": str(acquisition_datetime)
                    }
                elif ms_level == 1:
                    row_data = {
                        "SpectrumID": spectrum.getNativeID(),
                        "MSLevel": float(ms_level),
                        "Charge": None,
                        "MS_peaks": float(peak_per_ms),
                        "Base_Peak_Intensity": float(base_peak_intensity) if base_peak_intensity is not None else None,
                        "Summed_Peak_Intensities": float(total_intensity) if total_intensity is not None else None,
                        "Retention_Time": float(rt),
                        "Exp_Mass_To_Charge": None,
                        "AcquisitionDateTime": str(acquisition_datetime)
                    }
                else:
                    continue

                batch_data.append(row_data)

                # Write batch when it reaches specified size
                if len(batch_data) >= batch_size:
                    batch_table = pa.Table.from_pylist(batch_data, schema=parquet_schema)

                    if parquet_writer is None:
                        parquet_writer = pq.ParquetWriter(where=output_path, schema=parquet_schema, compression='gzip')

                    parquet_writer.write_table(batch_table)
                    batch_data = []

            # Write any remaining data
            if batch_data:
                batch_table = pa.Table.from_pylist(batch_data, schema=parquet_schema)
                if parquet_writer is None:
                    parquet_writer = pq.ParquetWriter(where=output_path, schema=parquet_schema, compression='gzip')
                parquet_writer.write_table(batch_table)

            # Write spectrum data if id_only
            if id_only and psm_parts:
                spectrum_table = pa.Table.from_pylist(
                    psm_parts,
                    schema=pa.schema([
                        ('scan', pa.string()),
                        ('ms_level', pa.int32()),
                        ('mz', pa.list_(pa.float64())),
                        ('intensity', pa.list_(pa.float64()))
                    ])
                )
                pq.write_table(
                    spectrum_table,
                    f"{Path(file_name).stem}_spectrum_df.parquet",
                    compression='gzip'
                )

            if parquet_writer is not None:
                parquet_writer.close()
            return output_path

        finally:
            if parquet_writer:
                parquet_writer.close()

    def batch_write_bruker_d(
        file_name: str,
        output_path: str,
        parquet_schema: pa.Schema,
        batch_size: int = 10000
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

            # Set up parquet writer
            parquet_writer = pq.ParquetWriter(
                output_path,
                parquet_schema,
                compression='gzip'
            )

            try:
                # Stream data in batches
                for chunk in pd.read_sql_query(
                    """
                    SELECT 
                        Id, 
                        CASE 
                            WHEN MsMsType IN (8, 9) THEN 2 
                            WHEN MsMsType = 0 THEN 1 
                            ELSE NULL 
                        END as MsMsType,
                        NumPeaks, 
                        MaxIntensity, 
                        SummedIntensities, 
                        Time,
                        Charge,
                        MonoisotopicMz
                    FROM frames
                    """,
                    conn,
                    chunksize=batch_size
                ):
                    # Prepare batch
                    chunk['AcquisitionDateTime'] = acquisition_date_time

                    # Convert to pyarrow table and write
                    batch_table = pa.Table.from_pandas(chunk, schema=parquet_schema)
                    parquet_writer.write_table(batch_table)

            finally:
                parquet_writer.close()

        return output_path

    # Resolve file path
    ms_path = _resolve_ms_path(ms_path)
    output_path = f"{Path(ms_path).stem}_ms_info.parquet"

    # Choose processing method based on file type
    if Path(ms_path).suffix == ".d":
        batch_write_bruker_d(file_name=ms_path, output_path=output_path, parquet_schema=schema, batch_size=batch_size)
    elif Path(ms_path).suffix.lower() in [".mzml"]:
        batch_write_mzml(file_name=ms_path, parquet_schema=schema, output_path=output_path, id_only=id_only, batch_size=batch_size)
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
    valid_extensions = {'.d', '.mzml', '.mzML'}
    candidates = [
        str(c.resolve())
        for c in candidates
        if c.suffix.lower() in valid_extensions
    ]

    if len(candidates) == 1:
        return candidates[0]

    raise FileNotFoundError(f"No unique file found for {ms_path}")

if __name__ == "__main__":
    mzml_statistics()