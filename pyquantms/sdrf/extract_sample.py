from pathlib import Path
import pandas as pd
import click

@click.command("extract_sample")
@click.option("--expdesign", help = "Experiment design file", type=click.Path(exists=True))
@click.pass_context
def extract_sample_from_expdesign(cxt, expdesign: str) -> None:
    """
    Extract sample information from an experiment design file
    :param cxt: Context object
    :param expdesign: Experiment design file
    :return: None
    """
    data = pd.read_csv(expdesign, sep="\t", header=0, dtype=str)
    fTable = data.dropna()

    # two table formats
    with open(expdesign, "r") as f:
        lines = f.readlines()
        empty_row = lines.index("\n")
        s_table = [i.replace("\n", "").split("\t") for i in lines[empty_row + 1:]][1:]
        s_header = lines[empty_row + 1].replace("\n", "").split("\t")
        s_DataFrame = pd.DataFrame(s_table, columns=s_header)

    sample_dt = pd.DataFrame()
    if "MSstats_Mixture" not in s_DataFrame.columns:
        fTable = fTable[["Spectra_Filepath", "Sample"]]
        fTable.to_csv(f"{Path(expdesign).stem}_sample.csv", sep="\t", index=False)
    else:
        fTable.drop_duplicates(subset=["Spectra_Filepath"], inplace=True)
        for _, row in fTable.iterrows():
            mixture_id = s_DataFrame[s_DataFrame["Sample"] == row["Sample"]]["MSstats_Mixture"]
            sample_dt = sample_dt.append({"Spectra_Filepath": row["Spectra_Filepath"], "Sample": mixture_id},
                                         ignore_index=True)
        sample_dt.to_csv(f"{Path(expdesign).stem}_sample.csv", sep="\t", index=False)