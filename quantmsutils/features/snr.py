import re

import click
import numpy as np
from pyopenms import MzMLFile, MSExperiment, SpectrumLookup
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
    "--idxml",
    help="The idxml file with the PSMs corresponding to the mzML file",
    required=True,
)
@click.option(
    "--output",
    help="The output idXML file with the signal-to-noise ratio",
    required=True,
)
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

    lookup = SpectrumLookup()
    lookup.readSpectra(mzml_file, "scan=(?<SCAN>\\d+)")

    protein_ids = []
    peptide_ids = []
    oms.IdXMLFile().load(idxml, protein_ids, peptide_ids)

    for peptide in peptide_ids:
        spectrum_reference = peptide.getMetaValue("spectrum_reference")
        scan_number = int(re.findall(r"(spectrum|scan)=(\d+)", spectrum_reference)[0][1])

        try:
           index = lookup.findByScanNumber(scan_number)
           spectrum = mzml_file.getSpectrum(index)
           intensity_array = spectrum.get_peaks()[1]
           snr = compute_signal_to_noise(intensity_array)
           for hit in peptide.getHits():
               hit.setMetaValue("quantms:SNR", str(snr))
        except IndexError:
            message = "scan_number" + str(scan_number) + "not found in file: " + ms_path
            print(message)

    # Add quantms:SNR as a feature
    search_parameters = protein_ids[0].getSearchParameters()
    features = search_parameters.getMetaValue("extra_features")
    extra_features = features + ",quantms:SNR"
    search_parameters.setMetaValue("extra_features", extra_features)
    protein_ids[0].setSearchParameters(search_parameters)

    oms.IdXMLFile().store(output, protein_ids, peptide_ids)
    print("The output file with the signal-to-noise ratio has been saved at: ", output)
