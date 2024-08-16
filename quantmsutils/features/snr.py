import re

import click
import numpy as np
from pyopenms import MzMLFile, MSExperiment, SpectrumLookup
import pyopenms as oms
from scipy.stats import entropy


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


@click.command("spectrum2feature")
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
def spectrum2feature(ctx, ms_path: str, idxml: str, output: str) -> None:
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

    result_peptides = []
    for peptide in peptide_ids:
        spectrum_reference = peptide.getMetaValue("spectrum_reference")
        scan_number = int(
            re.findall(r"(spectrum|scan)=(\d+)", spectrum_reference)[0][1]
        )

        try:
            index = lookup.findByScanNumber(scan_number)
            spectrum = mzml_file.getSpectrum(index)
            intensity_array = spectrum.get_peaks()[1]
            mz_array = spectrum.get_peaks()[0]

            ## Compute signal to noise ratio
            snr = compute_signal_to_noise(intensity_array)

            # Spectral Entropy
            tic = np.sum(intensity_array)
            normalized_intensities = intensity_array / tic
            spectral_entropy = entropy(normalized_intensities)

            # Fraction of TIC in Top 10 Peaks
            top_n_peaks = np.sort(intensity_array)[-10:]
            fraction_tic_top_10 = np.sum(top_n_peaks) / tic

            # Intensity Weighted m/z Standard Deviation
            weighted_mean_mz = np.sum(np.array(mz_array) * normalized_intensities)
            weighted_std_mz = np.sqrt(
                np.sum(
                    normalized_intensities
                    * (np.array(mz_array) - weighted_mean_mz) ** 2
                )
            )

            for hit in peptide.getHits():
                hit.setMetaValue("quantms:SNR", str(snr))
                hit.setMetaValue("quantms:SpectralEntropy", str(spectral_entropy))
                hit.setMetaValue(
                    "quantms:FracTICinTop10Peaks", str(fraction_tic_top_10)
                )
                hit.setMetaValue("quantms:WeightedStdMz", str(weighted_std_mz))

                # update hit in peptidehits list
                peptide.setHits([hit])
            ## update peptide in peptide_ids list
            result_peptides.append(peptide)

        except IndexError:
            message = "scan_number" + str(scan_number) + "not found in file: " + ms_path
            print(message)

    # Add quantms:SNR as a feature
    search_parameters = protein_ids[0].getSearchParameters()
    features = search_parameters.getMetaValue("extra_features")
    extra_features = (
        features
        + ",quantms:SNR,quantms:SpectralEntropy,quantms:FracTICinTop10Peaks,quantms:WeightedStdMz"
    )
    search_parameters.setMetaValue("extra_features", extra_features)
    protein_ids[0].setSearchParameters(search_parameters)

    oms.IdXMLFile().store(output, protein_ids, peptide_ids)
    print("The output file with the signal-to-noise ratio has been saved at: ", output)
