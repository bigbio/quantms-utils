import click
import numpy as np
from pyopenms import MzMLFile, MSExperiment
import pyopenms as oms


def compute_signal_to_noise(intensities):
    """
    Compute the signal-to-noise ratio for a given spectrum
    :param intensities: intensity values
    :return:
    """
    rmsd = np.sqrt(np.mean(np.square(intensities)))

    # Calculate SNR
    snr = np.max(intensities) / rmsd

    return snr

@click.command("snr")
@click.option("--ms_path", type=click.Path(exists=True), required=True)
@click.option(
    "--idxml", help="The idxml file with the PSMs corresponding to the mzML file", required=True
)
@click.option("--output", help="The output idXML file with the signal-to-noise ratio", required=True)
@click.pass_context
def snr_cmd(ctx, ms_path: str, idxml: str, output: str) -> None:
    """
    The snr function computes the signal-to-noise ratio for a given mass spectrometry file and its corresponding
    identification file.

    :param ctx: Click context
    :param ms_path: A string specifying the path to the mass spectrometry file.
    :param idxml: A string specifying the path to the identification file.
    :param output: A string specifying the path to the output file.
    """

    mzml_file = MSExperiment()
    MzMLFile().load(ms_path, mzml_file)

    protein_ids = []
    peptide_ids = []
    oms.IdXMLFile().load(idxml, protein_ids, peptide_ids)



    oms.IdXMLFile().store(output, protein_ids, peptide_ids)
    print("The output file with the signal-to-noise ratio has been saved at: ", output)



