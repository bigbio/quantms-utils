import re
import sqlite3
from pathlib import Path

import click
import pandas as pd
import pyarrow
from pyopenms import MSExperiment, MzMLFile


@click.command("mzmlstats")
@click.option("--ms_path", type=click.Path(exists=True))
@click.option(
    "--id_only", is_flag=True, help="Generate a csv with the spectrum id and the peaks"
)
@click.pass_context
def mzml_statistics(ctx, ms_path: str, id_only: bool = False) -> None:
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
    file_columns = [
        "SpectrumID",
        "MSLevel",
        "Charge",
        "MS_peaks",
        "Base_Peak_Intensity",
        "Summed_Peak_Intensities",
        "Retention_Time",
        "Exp_Mass_To_Charge",
        "AcquisitionDateTime",
    ]

    def parse_mzml(file_name: str, file_columns: list, id_only: bool = False):
        """
        Parse mzML file and return a pandas DataFrame with the information. If id_only is True, it will also save a csv.
        @param file_name: The file name of the mzML file
        @param file_columns: The columns of the DataFrame
        @param id_only: If True, it will save a csv with the spectrum id, mz and intensity
        @return: A pandas DataFrame with the information of the mzML file
        """

        info = []
        psm_part_info = []
        exp = MSExperiment()
        acquisition_datetime = exp.getDateTime().get()
        MzMLFile().load(file_name, exp)
        for spectrum in exp:
            id_ = spectrum.getNativeID()
            ms_level = spectrum.getMSLevel()
            rt = spectrum.getRT() if spectrum.getRT() else None

            peaks_tuple = spectrum.get_peaks()
            peak_per_ms = len(peaks_tuple[0])

            if not spectrum.metaValueExists("base peak intensity"):
                bpc = max(peaks_tuple[1]) if len(peaks_tuple[1]) > 0 else None
            else:
                bpc = spectrum.getMetaValue("base peak intensity")

            if not spectrum.metaValueExists("total ion current"):
                tic = sum(peaks_tuple[1]) if len(peaks_tuple[1]) > 0 else None
            else:
                tic = spectrum.getMetaValue("total ion current")

            if ms_level == 1:
                info_list = [
                    id_,
                    ms_level,
                    None,
                    peak_per_ms,
                    bpc,
                    tic,
                    rt,
                    None,
                    acquisition_datetime,
                ]
            elif ms_level == 2:
                charge_state = spectrum.getPrecursors()[0].getCharge()
                emz = (
                    spectrum.getPrecursors()[0].getMZ()
                    if spectrum.getPrecursors()[0].getMZ()
                    else None
                )
                info_list = [
                    id_,
                    ms_level,
                    charge_state,
                    peak_per_ms,
                    bpc,
                    tic,
                    rt,
                    emz,
                    acquisition_datetime,
                ]
                mz_array = peaks_tuple[0]
                intensity_array = peaks_tuple[1]
            else:
                info_list = [
                    id_,
                    ms_level,
                    None,
                    None,
                    None,
                    None,
                    rt,
                    None,
                    acquisition_datetime,
                ]

            if id_only and ms_level == 2:
                psm_part_info.append(
                    [
                        re.findall(r"[scan|spectrum]=(\d+)", id_)[0],
                        ms_level,
                        mz_array,
                        intensity_array,
                    ]
                )
            info.append(info_list)

        if id_only and len(psm_part_info) > 0:
            pd.DataFrame(
                psm_part_info, columns=["scan", "ms_level", "mz", "intensity"]
            ).to_parquet(
                f"{Path(ms_path).stem}_spectrum_df.parquet",
                index=False,
                compression="gzip",
            )

        return pd.DataFrame(info, columns=file_columns)

    def parse_bruker_d(file_name: str, file_columns: list):
        sql_filepath = f"{file_name}/analysis.tdf"
        if not Path(sql_filepath).exists():
            msg = f"File '{sql_filepath}' not found"
            raise FileNotFoundError(msg)
        conn = sqlite3.connect(sql_filepath)
        c = conn.cursor()

        datetime_cmd = (
            "SELECT Value FROM GlobalMetadata WHERE key='AcquisitionDateTime'"
        )
        acquisition_date_time = c.execute(datetime_cmd).fetchall()[0][0]

        df = pd.read_sql_query(
            "SELECT Id, MsMsType, NumPeaks, MaxIntensity, SummedIntensities, Time FROM frames",
            conn,
        )
        df["AcquisitionDateTime"] = acquisition_date_time

        # {8:'DDA-PASEF', 9:'DIA-PASEF'}
        if 8 in df["MsMsType"].values:
            mslevel_map = {0: 1, 8: 2}
        elif 9 in df["MsMsType"].values:
            mslevel_map = {0: 1, 9: 2}
        else:
            msg = f"Unrecognized ms type '{df['MsMsType'].values}'"
            raise ValueError(msg)
        df["MsMsType"] = df["MsMsType"].map(mslevel_map)

        try:
            # This line raises an sqlite error if the table does not exist
            _ = conn.execute("SELECT * from Precursors LIMIT 1").fetchall()
            precursor_df = pd.read_sql_query("SELECT * from Precursors", conn)
        except sqlite3.OperationalError as e:
            if "no such table: Precursors" in str(e):
                print(
                    f"No precursors recorded in {file_name}, This is normal for DIA data."
                )
                precursor_df = pd.DataFrame()
            else:
                raise

        if len(df) == len(precursor_df):
            df = pd.concat([df, precursor_df["Charge", "MonoisotopicMz"]], axis=1)
            df["Charge"] = df["Charge"].fillna(0)
        else:
            df[["Charge", "Exp_Mass_To_Charge"]] = None, None

        df = df[
            [
                "Id",
                "MsMsType",
                "Charge",
                "NumPeaks",
                "MaxIntensity",
                "SummedIntensities",
                "Time",
                "Exp_Mass_To_Charge",
                "AcquisitionDateTime",
            ]
        ]
        df.columns = pd.Index(file_columns)

        return df

    if not (Path(ms_path).exists()):
        print(f"Not found '{ms_path}', trying to find alias")
        ms_path_path = Path(ms_path)
        path_stem = str(ms_path_path.stem)
        candidates = (
            list(ms_path_path.parent.glob("*.d"))
            + list(ms_path_path.parent.glob("*.mzml"))
            + list(ms_path_path.parent.glob("*.mzML"))
        )

        candidates = [c for c in candidates if path_stem in str(c)]

        if len(candidates) == 1:
            ms_path = str(candidates[0].resolve())
        else:
            raise FileNotFoundError()

    if Path(ms_path).suffix == ".d" and Path(ms_path).is_dir():
        ms_df = parse_bruker_d(ms_path, file_columns)
    elif Path(ms_path).suffix in [".mzML", ".mzml"]:
        ms_df = parse_mzml(ms_path, file_columns, id_only)
    else:
        msg = f"Unrecognized or the mass spec file '{ms_path}' do not exist"
        raise RuntimeError(msg)

    ms_df.to_parquet(
        f"{Path(ms_path).stem}_ms_info.parquet",
        engine="pyarrow",
        index=False,
        compression="gzip",
    )
