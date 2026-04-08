"""
Convert DIA-NN output to MSstats format.
License: Apache 2.0
Authors: Hong Wong, Yasset Perez-Riverol
Revisions:
    2023-Aug-05: J. Sebastian Paez
"""

import logging
from pathlib import Path

import click
import pandas as pd
import pyarrow.parquet as pq
from pyopenms import AASequence

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
REVISION = "0.1.1"

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


@click.command("diann2msstats", short_help="Convert DIA-NN output to MSstats")
@click.option("--report", "-r", "report_path", required=True, type=click.Path(exists=True))
@click.option("--exp_design", "-d", required=True, type=click.Path(exists=True))
@click.option("--qvalue_threshold", "-q", type=float, required=True)
def diann2msstats(
    report_path,
    exp_design,
    qvalue_threshold,
):
    """
    Convert DIA-NN output to MSstats format for downstream analysis.

    :param report_path: DIA-NN main report file (.tsv or .parquet)
    :param exp_design: Experimental design file
    :param qvalue_threshold: Q-value filter threshold
    """
    logger.debug(f"Revision {REVISION}")
    logger.debug("Reading input files...")
    report = load_report(report_path, qvalue_threshold)
    s_data_frame, f_table = get_exp_design_dfs(exp_design)

    msstats_columns_keep = [
        "Protein.Names",
        "Modified.Sequence",
        "Precursor.Charge",
        "Precursor.Quantity",
        "Run",
    ]

    out_msstats_columns = [
        "ProteinName",
        "PeptideSequence",
        "PrecursorCharge",
        "Intensity",
        "Run",
    ]

    if "Channel" in report.columns and report["Channel"].nunique() > 1:
        logger.info("Multiplexing detected. Applying label mapping...")
        msstats_columns_keep.append("Channel")
        out_msstats_columns.append("IsotopeLabelType")

    logger.debug("Converting to MSstats format...")
    if "Decoy" in report.columns:
        out_msstats = report[report["Decoy"] != 1][msstats_columns_keep].copy()
    else:
        out_msstats = report[msstats_columns_keep].copy()

    out_msstats.columns = out_msstats_columns
    out_msstats = out_msstats[out_msstats["Intensity"] != 0]

    out_msstats["PeptideSequence"] = out_msstats["PeptideSequence"].apply(_sanitize_sequence)
    out_msstats.loc[:, "PeptideSequence"] = out_msstats.apply(
        lambda x: (
            AASequence.fromString(x["PeptideSequence"]).toString()
            if "^" not in x["PeptideSequence"]
            else "^" + AASequence.fromString(x["PeptideSequence"].replace("^", "")).toString()
        ),
        axis=1,
    )
    out_msstats["FragmentIon"] = "NA"
    out_msstats["ProductCharge"] = "0"

    if "IsotopeLabelType" in out_msstats.columns:
        out_msstats = out_msstats[out_msstats["IsotopeLabelType"].notna()]
        out_msstats = out_msstats[out_msstats["IsotopeLabelType"].str.strip() != ""]

        # is multiplexing
        f_table_cols = ["Fraction", "Sample", "run", "Label"]
        merge_keys = ["Run", "IsotopeLabelType"]
    else:
        out_msstats["IsotopeLabelType"] = "L"

        f_table_cols = ["Fraction", "Sample", "run"]
        merge_keys = ["Run"]

    logger.debug(f"out_msstats ({out_msstats.shape}) >>>")
    logger.debug("Adding Fraction, BioReplicate, Condition columns")

    design_lookup = (
        s_data_frame[["Sample", "MSstats_Condition", "MSstats_BioReplicate"]]
        .merge(f_table[f_table_cols], on="Sample")
        .rename(
            columns={
                "run": "Run",
                "MSstats_BioReplicate": "BioReplicate",
                "MSstats_Condition": "Condition",
                "Label": "IsotopeLabelType",
            },
            errors="ignore"
        )
        .drop(columns=["Sample"])
    )
    out_msstats = out_msstats.merge(design_lookup, on=merge_keys, how="left", validate="many_to_one")

    unmatched = out_msstats["BioReplicate"].isna()
    if unmatched.any():
        bad_runs = out_msstats.loc[unmatched, "Run"].unique().tolist()
        logger.warning(
            "Run(s) in DIA-NN report have no match in experimental design: %s. "
            "These rows will be dropped. Check that Run names (spectra file stems) match Spectra_Filepath in the design.",
            bad_runs,
        )
        out_msstats = out_msstats.dropna(subset=["BioReplicate"])
    exp_out_prefix = Path(exp_design).stem
    out_msstats.to_csv(exp_out_prefix + "_msstats_in.csv", sep=",", index=False)
    logger.info(f"MSstats input file is saved as {exp_out_prefix}_msstats_in.csv")


def _true_stem(x):
    """Return the file name stem (without extension)."""
    if x.endswith(".d.zip"):
        return Path(Path(x).stem).stem
    else:
        return Path(x).stem


def get_exp_design_dfs(exp_design_file):
    logger.info(f"Reading experimental design file: {exp_design_file}")
    with open(exp_design_file, "r") as f:
        data = [line.replace("\r\n", "\n").replace("\r", "\n") for line in f.readlines()]

    # Auto-detect format: new unified flat table (from convert-diann) vs legacy two-table
    # Prefer unified parsing whenever the header matches unified columns,
    # regardless of extra blank/whitespace-only lines elsewhere in the file.
    header_line = data[0].replace("\n", "") if data else ""
    is_unified_format = "Condition" in header_line and "BioReplicate" in header_line and "Filename" in header_line

    if is_unified_format:
        return _parse_unified_design(exp_design_file)
    else:
        return _parse_legacy_design(data, exp_design_file)


def _parse_unified_design(exp_design_file):
    """Parse the unified flat diann_design.tsv format from convert-diann.

    Returns (s_data_frame, f_table) matching the legacy interface:
    - s_data_frame: DataFrame with Sample, MSstats_Condition, MSstats_BioReplicate
    - f_table: DataFrame with Fraction, Sample, run (+ other columns)
    """
    df = pd.read_csv(exp_design_file, sep="\t")

    # Validate required columns
    required_cols = {
        "Filename", "Fraction", "Sample", "Condition", "BioReplicate", "Label", "LabelType"
    }

    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Unified design file is missing required columns: {', '.join(sorted(missing))}. "
            f"Expected: {', '.join(sorted(required_cols))}"
        )

    df["run"] = df["Filename"].apply(_true_stem)

    # Multiplexing
    if df["Label"].nunique() > 1:

        if "silac" in df["LabelType"].values:
            silac_dict = {
                "SILAC light": "L",
                "SILAC medium": "M",
                "SILAC heavy": "H",
            }
            df["Label"] = df["Label"].replace(silac_dict)

        f_table = df[["Filename", "Fraction", "Sample", "run", "Label"]].copy()
    else:
        f_table = df[["Filename", "Fraction", "Sample", "run"]].copy()

    # Validate each Sample maps to exactly one (Condition, BioReplicate) pair
    unique_mapping = df[["Sample", "Condition", "BioReplicate"]].drop_duplicates()
    duplicated_samples = unique_mapping["Sample"][unique_mapping["Sample"].duplicated(keep=False)].unique()
    if len(duplicated_samples) > 0:
        raise ValueError(
            "Inconsistent experimental design: Sample(s) "
            f"{', '.join(map(str, duplicated_samples))} map to multiple "
            "(Condition, BioReplicate) combinations."
        )

    s_data_frame = unique_mapping.rename(
        columns={"Condition": "MSstats_Condition", "BioReplicate": "MSstats_BioReplicate"}
    )

    return s_data_frame, f_table


def _parse_legacy_design(data, exp_design_file):
    """Parse the legacy two-table experimental design format (blank line separator)."""
    try:
        empty_row = data.index("\n")
    except ValueError:
        raise ValueError(
            f"Could not find blank separator row in {exp_design_file}. "
            "Ensure the file contains a blank line between the file and sample tables."
        ) from None
    f_table = [i.replace("\n", "").split("\t") for i in data[1:empty_row]]
    f_header = data[0].replace("\n", "").split("\t")
    f_table = pd.DataFrame(f_table, columns=f_header)
    f_table.loc[:, "run"] = f_table.apply(lambda x: _true_stem(x["Spectra_Filepath"]), axis=1)

    s_table = [i.replace("\n", "").split("\t") for i in data[empty_row + 1 :]][1:]
    s_header = data[empty_row + 1].replace("\n", "").split("\t")
    s_data_frame = pd.DataFrame(s_table, columns=s_header)

    return s_data_frame, f_table


def load_report(report_path, qvalue_threshold: float) -> pd.DataFrame:
    """Load DIA-NN report from TSV or Parquet, detecting format from file extension."""
    path = Path(report_path)
    remain_cols = [
        "Run",
        "Protein.Names",
        "Modified.Sequence",
        "Precursor.Charge",
        "Precursor.Quantity",
        "Q.Value",
    ]
    if path.suffix == ".parquet":
        pq_schema = pq.read_schema(path)
        use_cols = remain_cols + [col for col in ["Decoy", "Channel"] if col in pq_schema.names]
        report = pd.read_parquet(path, columns=use_cols)
    else:
        tsv_header = pd.read_csv(path, sep="\t", header=0, nrows=0).columns.tolist()
        tsv_cols = remain_cols + (["Decoy"] if "Decoy" in tsv_header else [])
        report = pd.read_csv(path, sep="\t", header=0, usecols=tsv_cols)

    report = report[report["Q.Value"] < qvalue_threshold]
    return report


def _sanitize_sequence(seq):
    seq = seq.replace("(SILAC)", "")
    return seq
