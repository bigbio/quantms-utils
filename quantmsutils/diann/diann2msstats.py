"""
Convert DIA-NN output to MSstats format.
License: Apache 2.0
Authors: Hong Wong, Yasset Perez-Riverol
Revisions:
    2023-Aug-05: J. Sebastian Paez
"""

import logging
import os
from pathlib import Path
from typing import List

import click
import pandas as pd
from pyopenms import AASequence
from pyopenms.Constants import PROTON_MASS_U

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)
pd.set_option("display.width", 1000)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
REVISION = "0.1.1"

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


@click.command("diann2msstats", short_help="Convert DIA-NN output to MSstats")
@click.option("--folder", "-f")
@click.option("--exp_design", "-d")
@click.option("--diann_version", "-v")
@click.option("--dia_params", "-p")
@click.option("--charge", "-c")
@click.option("--missed_cleavages", "-m")
@click.option("--qvalue_threshold", "-q", type=float)
def diann2msstats(
    folder,
    exp_design,
    dia_params,
    diann_version,
    charge,
    missed_cleavages,
    qvalue_threshold,
):
    """
    Convert DIA-NN output to MSstats format for downstream analysis.

    :param folder: Folder containing DIA-NN main report, experimental design,
        protein sequence FASTA, version file, and ms_info parquet
    :param exp_design: Experimental design file
    :param dia_params: DIA parameters (semicolon-separated)
    :param diann_version: Path to DIA-NN version file
    :param charge: Max precursor charge from DIA-NN
    :param missed_cleavages: Allowed missed cleavages
    :param qvalue_threshold: Q-value filter threshold
    """
    logger.debug(f"Revision {REVISION}")
    logger.debug("Reading input files...")
    diann_directory = DiannDirectory(folder, diann_version_file=diann_version)
    report = diann_directory.main_report_df(qvalue_threshold=qvalue_threshold)
    s_data_frame, f_table = get_exp_design_dfs(exp_design)

    msstats_columns_keep = [
        "Protein.Names",
        "Modified.Sequence",
        "Precursor.Charge",
        "Precursor.Quantity",
        "Run",
    ]

    logger.debug("Converting to MSstats format...")
    if "Decoy" in report.columns:
        out_msstats = report[report["Decoy"] != 1][msstats_columns_keep]
    else:
        out_msstats = report[msstats_columns_keep]

    out_msstats.columns = [
        "ProteinName",
        "PeptideSequence",
        "PrecursorCharge",
        "Intensity",
        "Run",
    ]
    out_msstats = out_msstats[out_msstats["Intensity"] != 0]

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
    out_msstats["IsotopeLabelType"] = "L"

    logger.debug(f"out_msstats ({out_msstats.shape}) >>>")
    logger.debug("Adding Fraction, BioReplicate, Condition columns")

    out_msstats = out_msstats.merge(
        (
            s_data_frame[["Sample", "MSstats_Condition", "MSstats_BioReplicate"]]
            .merge(f_table[["Fraction", "Sample", "run"]], on="Sample")
            .rename(
                columns={
                    "run": "Run",
                    "MSstats_BioReplicate": "BioReplicate",
                    "MSstats_Condition": "Condition",
                }
            )
            .drop(columns=["Sample"])
        ),
        on="Run",
        validate="many_to_one",
    )
    exp_out_prefix = Path(exp_design).stem
    out_msstats.to_csv(exp_out_prefix + "_msstats_in.csv", sep=",", index=False)
    logger.info(f"MSstats input file is saved as {exp_out_prefix}_msstats_in.csv")


def _true_stem(x):
    """Return the file name stem (without extension)."""
    return os.path.basename(x).split(".")[0]


def get_exp_design_dfs(exp_design_file):
    logger.info(f"Reading experimental design file: {exp_design_file}")
    with open(exp_design_file, "r") as f:
        data = f.readlines()
        empty_row = data.index("\n")
        f_table = [i.replace("\n", "").split("\t") for i in data[1:empty_row]]
        f_header = data[0].replace("\n", "").split("\t")
        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table.loc[:, "run"] = f_table.apply(lambda x: _true_stem(x["Spectra_Filepath"]), axis=1)

        s_table = [i.replace("\n", "").split("\t") for i in data[empty_row + 1 :]][1:]
        s_header = data[empty_row + 1].replace("\n", "").split("\t")
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    return s_data_frame, f_table


def compute_mass_modified_peptide(peptide_seq: str) -> float:
    """Compute peptide mass including modifications using pyopenms."""
    peptide_parts: List[str] = []
    not_mod = True
    aa_mass = {
        "X": "X[178.98493453312]",
        "U": "X[132.94306553312]",
        "O": "X[237.14773053312]",
    }
    valid_aa = set("GAVLIFMPWSTYNQDEKRH")
    for aa in peptide_seq:
        if aa == "^":
            continue
        if aa == "(":
            not_mod = False
        elif aa == ")":
            not_mod = True
        if aa in aa_mass and not_mod:
            aa = aa_mass[aa]
        elif aa not in valid_aa and not_mod and aa != ")":
            logger.info(f"Unknown amino acid with mass not known:{aa}")
        peptide_parts.append(aa)
    new_peptide_seq = "".join(peptide_parts)
    return AASequence.fromString(new_peptide_seq).getMonoWeight()


class DiannDirectory:
    def __init__(self, base_path, diann_version_file):
        self.base_path = Path(base_path)
        if not self.base_path.exists() or not self.base_path.is_dir():
            raise NotADirectoryError(f"Path {self.base_path} does not exist or is not a directory")
        self.diann_version_file = Path(diann_version_file)
        if not self.diann_version_file.is_file():
            raise FileNotFoundError(f"Path {self.diann_version_file} does not exist")

    def find_first_file_with_suffix(self, suffix: str) -> os.PathLike:
        try:
            return next(self.base_path.glob(f"**/*{suffix}"))
        except StopIteration:
            raise FileNotFoundError(f"Could not find file with suffix {suffix}")

    @property
    def report(self) -> os.PathLike:
        try:
            return self.find_first_file_with_suffix("report.tsv")
        except FileNotFoundError:
            return self.find_first_file_with_suffix("report.parquet")

    @property
    def diann_version(self) -> str:
        logger.debug("Validating DIANN version")
        diann_version_id = None
        with open(self.diann_version_file) as f:
            for line in f:
                if "DIA-NN" in line:
                    diann_version_id = line.rstrip("\n").split(": ")[1]
                    break
        if diann_version_id is None:
            raise ValueError(f"Could not find DIA-NN version in file {self.diann_version_file}")
        return diann_version_id

    def main_report_df(self, qvalue_threshold: float) -> pd.DataFrame:
        remain_cols = [
            "Run",
            "Protein.Group",
            "Protein.Names",
            "Protein.Ids",
            "PG.MaxLFQ",
            "RT",
            "Global.Q.Value",
            "Lib.Q.Value",
            "PEP",
            "Precursor.Normalised",
            "Precursor.Id",
            "Q.Value",
            "Modified.Sequence",
            "Stripped.Sequence",
            "Precursor.Charge",
            "Precursor.Quantity",
            "Global.PG.Q.Value",
        ]
        if self.diann_version.startswith("1."):
            report = pd.read_csv(
                self.report, sep="\t", header=0, usecols=remain_cols + ["MS2.Scan"]
            )
        else:
            report = pd.read_parquet(self.report, columns=remain_cols + ["Decoy"])

        report = report[report["Q.Value"] < qvalue_threshold]

        uniq_masses = {
            k: compute_mass_modified_peptide(k) for k in report["Modified.Sequence"].unique()
        }
        mass_vector = report["Modified.Sequence"].map(uniq_masses)
        report["Calculate.Precursor.Mz"] = (
            mass_vector + (PROTON_MASS_U * report["Precursor.Charge"])
        ) / report["Precursor.Charge"]

        precursor_index_map = {k: i for i, k in enumerate(report["Precursor.Id"].unique())}
        report["precursor.Index"] = report["Precursor.Id"].map(precursor_index_map)

        return report
