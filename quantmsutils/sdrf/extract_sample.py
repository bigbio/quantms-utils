import logging
from pathlib import Path

import click
import pandas as pd

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


@click.command(
    "openms2sample",
    short_help="Extract sample information from an experiment design file",
)
@click.option("--expdesign", help="Experiment design file", type=click.Path(exists=True))
@click.pass_context
def extract_sample_from_expdesign(cxt, expdesign: str) -> None:
    """
    Extract sample information from an experiment design file
    :param cxt: Context object
    :param expdesign: Experiment design file
    :return: None
    """
    data = pd.read_csv(expdesign, sep="\t", header=0, dtype=str)
    f_table = data.dropna()

    # two table formats
    with open(expdesign, "r") as f:
        lines = f.readlines()
        empty_row = lines.index("\n")
        s_table = [i.replace("\n", "").split("\t") for i in lines[empty_row + 1 :]][1:]
        s_header = lines[empty_row + 1].replace("\n", "").split("\t")
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    sample_dt = pd.DataFrame()
    if "MSstats_Mixture" not in s_data_frame.columns:
        f_table = f_table[["Spectra_Filepath", "Sample"]]
        f_table.to_csv(f"{Path(expdesign).stem}_sample.csv", sep="\t", index=False)
    else:
        f_table.drop_duplicates(subset=["Spectra_Filepath"], inplace=True)
        for _, row in f_table.iterrows():
            mixture_id = s_data_frame[s_data_frame["Sample"] == row["Sample"]][
                "MSstats_Mixture"
            ].values[0]
            sample_dt = sample_dt.append(
                {"Spectra_Filepath": row["Spectra_Filepath"], "Sample": mixture_id},
                ignore_index=True,
            )
        sample_dt.to_csv(f"{Path(expdesign).stem}_sample.csv", sep="\t", index=False)
