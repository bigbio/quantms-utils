# Written by Jonas Scheid under the MIT license
# Contributions by Yasset Perez-Riverol and Dai Chengxin
# This script is part of the quantmsutils package

import importlib.resources
import json
import logging
from typing import List

import click
import pyopenms as oms
from ms2rescore import package_data, rescore
from psm_utils import PSMList
from psm_utils.io.idxml import IdXMLReader, IdXMLWriter

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def parse_cli_arguments_to_config(
    config_file: str = None,
    feature_generators: str = None,
    ms2pip_model_dir: str = None,
    ms2pip_model: str = None,
    ms2_tolerance: float = None,
    calibration_set_size: float = None,
    rescoring_engine: str = None,
    rng: int = None,
    test_fdr: float = None,
    processes: int = None,
    spectrum_path: str = None,
    fasta_file: str = None,
    id_decoy_pattern: str = None,
    lower_score_is_better: bool = None,
    output_path: str = None,
    log_level: str = None,
    spectrum_id_pattern: str = None,
    psm_id_pattern: str = None
) -> dict:
    if config_file is None:
        config = json.load(
            importlib.resources.open_text(package_data, "config_default.json")
        )
    else:
        with open(config_file) as f:
            config = json.load(f)
    if feature_generators is not None:
        feature_generators_list = feature_generators.split(",")
        config["ms2rescore"]["feature_generators"] = {}
        if "basic" in feature_generators_list:
            config["ms2rescore"]["feature_generators"]["basic"] = {}
        if "ms2pip" in feature_generators_list:
            config["ms2rescore"]["feature_generators"]["ms2pip"] = {
                "model_dir": ms2pip_model_dir,
                "model": ms2pip_model,
                "ms2_tolerance": ms2_tolerance,
            }
        if "deeplc" in feature_generators_list:
            config["ms2rescore"]["feature_generators"]["deeplc"] = {
                "deeplc_retrain": False,
                "calibration_set_size": calibration_set_size,
            }
        if "maxquant" in feature_generators_list:
            config["ms2rescore"]["feature_generators"]["maxquant"] = {}
        if "ionmob" in feature_generators:
            config["ms2rescore"]["feature_generators"]["ionmob"] = {}

    if rescoring_engine is not None:
        # Reset rescoring engine dict we want to allow only computing features
        config["ms2rescore"]["rescoring_engine"] = {}
        if rescoring_engine == "mokapot":
            config["ms2rescore"]["rescoring_engine"]["mokapot"] = {
                "write_weights": True,
                "write_txt": False,
                "write_flashlfq": False,
                "rng": rng,
                "test_fdr": test_fdr,
                "max_workers": processes,
            }
        if rescoring_engine == "percolator":
            logging.info(
                "Percolator rescoring engine has been specified. Use the idXML containing rescoring features and run Percolator in a separate step."
            )

    if ms2pip_model_dir is not None:
        config["ms2rescore"]["ms2pip_model_dir"] = ms2pip_model_dir
    if ms2pip_model is not None:
        config["ms2rescore"]["ms2pip_model"] = ms2pip_model
    if ms2_tolerance is not None:
        config["ms2rescore"]["ms2_tolerance"] = ms2_tolerance
    if calibration_set_size is not None:
        config["ms2rescore"]["calibration_set_size"] = calibration_set_size
    if rng is not None:
        config["ms2rescore"]["rng"] = rng
    if spectrum_path is not None:
        config["ms2rescore"]["spectrum_path"] = spectrum_path
    if fasta_file is not None:
        config["ms2rescore"]["fasta_file"] = fasta_file
    if id_decoy_pattern is not None:
        config["ms2rescore"]["id_decoy_pattern"] = id_decoy_pattern
    if lower_score_is_better is not None:
        config["ms2rescore"]["lower_score_is_better"] = lower_score_is_better
    if processes is None:
        processes = 1  # Default to single process
    config["ms2rescore"]["processes"] = processes
    if output_path is not None:
        config["ms2rescore"]["output_path"] = output_path
    else:
        raise ValueError("Output path must be specified.")
    if log_level is not None:
        config["ms2rescore"]["log_level"] = log_level
    if spectrum_id_pattern is not None:
        config["ms2rescore"]["spectrum_id_pattern"] = spectrum_id_pattern
    if psm_id_pattern is not None:
        config["ms2rescore"]["psm_id_pattern"] = psm_id_pattern

    return config


def rescore_idxml(input_file, output_file, config) -> None:
    """Rescore PSMs in an idXML file and keep other information unchanged."""
    # Read PSMs
    reader = IdXMLReader(input_file)
    psm_list = reader.read_file()

    # Rescore
    rescore(config, psm_list)

    # Filter out PeptideHits within PeptideIdentification(s) that could not be processed by all feature generators
    peptide_ids_filtered = filter_out_artifact_psms(psm_list, reader.peptide_ids)

    # Write
    writer = IdXMLWriter(output_file, reader.protein_ids, peptide_ids_filtered)
    writer.write_file(psm_list)


def filter_out_artifact_psms(
    psm_list: PSMList, peptide_ids: List[oms.PeptideIdentification]
) -> List[oms.PeptideIdentification]:
    """Filter out PeptideHits that could not be processed by all feature generators"""
    num_mandatory_features = max([len(psm.rescoring_features) for psm in psm_list])
    new_psm_list = PSMList(
        psm_list=[
            psm
            for psm in psm_list
            if len(psm.rescoring_features) == num_mandatory_features
        ]
    )

    # get differing peptidoforms of both psm lists
    psm_list_peptides = set(
        [next(iter(psm.provenance_data.items()))[1] for psm in psm_list]
    )
    new_psm_list_peptides = set(
        [next(iter(psm.provenance_data.items()))[1] for psm in new_psm_list]
    )
    not_supported_peptides = psm_list_peptides - new_psm_list_peptides

    # no need to filter if all peptides are supported
    if len(not_supported_peptides) == 0:
        return peptide_ids
    # Create new peptide ids and filter out not supported peptides
    new_peptide_ids = []
    for peptide_id in peptide_ids:
        new_hits = []
        for hit in peptide_id.getHits():
            if hit.getSequence().toString() in not_supported_peptides:
                continue
            new_hits.append(hit)
        if len(new_hits) == 0:
            continue
        peptide_id.setHits(new_hits)
        new_peptide_ids.append(peptide_id)
    logging.info(
        f"Removed {len(psm_list_peptides) - len(new_psm_list_peptides)} PSMs. Peptides not supported: {not_supported_peptides}"
    )
    return new_peptide_ids


@click.command(
    "ms2rescore",
    short_help="Rescore PSMs in an idXML file and keep other information unchanged.",
)
@click.option(
    "-p",
    "--psm_file",
    help="Path to PSM file (idXML)",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-s",
    "--spectrum_path",
    help="Path to MGF/mzML spectrum file or directory with spectrum files (default: derived from identification file)",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output_path",
    help="Path and stem for output file names (default: derive from identification file)",
)
@click.option(
    "-l", "--log_level", help="Logging level (default: `info`)", default="info"
)
@click.option(
    "-n",
    "--processes",
    help="Number of parallel processes available to MS²Rescore",
    type=int,
    default=16,
)
@click.option("-f", "--fasta_file", help="Path to FASTA file")
@click.option(
    "-t",
    "--test_fdr",
    help="The false-discovery rate threshold at which to evaluate the learned models. (default: 0.05)",
    default=0.05,
)
@click.option(
    "-fg",
    "--feature_generators",
    help="Comma-separated list of feature generators to use (default: `ms2pip,deeplc`). See rescoring doc for further information",
    default="",
)
@click.option(
    "-pipm",
    "--ms2pip_model",
    help="MS²PIP model (default: `Immuno-HCD`)",
    type=str,
    default="Immuno-HCD",
)
@click.option(
    "-md",
    "--ms2pip_model_dir",
    help="The path of MS²PIP model (default: `./`)",
    type=str,
    default="./",
)
@click.option(
    "-ms2tol",
    "--ms2_tolerance",
    help="Fragment mass tolerance [Da](default: `0.02`)",
    type=float,
    default=0.02,
)
@click.option(
    "-cs",
    "--calibration_set_size",
    help="Percentage of number of calibration set for DeepLC (default: `0.15`)",
    default=0.15,
)
@click.option(
    "-re",
    "--rescoring_engine",
    help="Either mokapot or percolator (default: `percolator`)",
    default="percolator",
    type=click.Choice(["mokapot", "percolator"]),
)
@click.option(
    "-rng",
    "--rng",
    help="Seed for mokapot's random number generator (default: `4711`)",
    type=int,
    default=4711,
)
@click.option(
    "-d",
    "--id_decoy_pattern",
    help="Regex decoy pattern (default: `DECOY_`)",
    default="^DECOY_",
)
@click.option(
    "-lsb",
    "--lower_score_is_better",
    help="Interpretation of primary search engine score (default: True)",
    default=True,
)
@click.option(
    "--config_file",
    help="Path to MS²Rescore config file (default: `config_default.json`)",
    default=None,
)
@click.option(
    "--spectrum_id_pattern",
    help="Regex pattern to extract index or scan number from spectrum file. Requires at least one capturing group.",
    default="(.*)",
)
@click.option(
    "--psm_id_pattern",
    help="Regex pattern to extract index or scan number from PSM file. Requires at least one capturing group.",
    default="(.*)",
)
@click.pass_context
def ms2rescore(
    ctx,
    psm_file: str,
    spectrum_path,
    output_path: str,
    log_level,
    processes,
    fasta_file,
    test_fdr,
    feature_generators,
    ms2pip_model_dir,
    ms2pip_model,
    ms2_tolerance,
    calibration_set_size,
    rescoring_engine,
    rng,
    id_decoy_pattern,
    lower_score_is_better,
    config_file: str,
    spectrum_id_pattern: str,
    psm_id_pattern: str
):
    """
    Rescore PSMs in an idXML file and keep other information unchanged.
    :param ms2pip_model_dir: Folder for models.
    :param ctx: Click context object
    :param psm_file: PSM file (idXML)
    :param spectrum_path: Spectrum file or dictionary with spectrum files (MGF/mzML)
    :param output_path: Output path for the new featured idXML file
    :param log_level: log_level for the logger
    :param processes: Number of parallel processes available to MS²Rescore
    :param fasta_file: Fasta file for the database search
    :param test_fdr: test FDR for the rescoring engine
    :param feature_generators: feature generators to use
    :param ms2pip_model: ms2pip model to use
    :param ms2_tolerance: ms2 tolerance
    :param calibration_set_size: calibration set size
    :param rescoring_engine: rescoring engine to use (mokapot or percolator)
    :param rng: random number generator seed
    :param id_decoy_pattern: id decoy pattern
    :param lower_score_is_better: lower score is better
    :param config_file: config file
    :param spectrum_id_pattern:egex pattern to extract index or scan number from spectrum file
    :param psm_id_pattern: Regex pattern to extract index or scan number from PSM file
    :return:
    """
    logging.getLogger().setLevel(log_level.upper())

    if output_path is None:
        output_path = psm_file.replace(".idXML", "_ms2rescore.idXML")

    if rescoring_engine == "moakapot":
        logging.warning(
            "Mokapot rescoring engine is not supported in this version. Please use Percolator."
        )
        raise ValueError(
            "Mokapot rescoring engine is not supported in this version. Please use Percolator."
        )

    config = parse_cli_arguments_to_config(
        config_file=config_file,
        output_path=output_path,
        feature_generators=feature_generators,
        ms2pip_model_dir=ms2pip_model_dir,
        ms2pip_model=ms2pip_model,
        processes=processes,
        ms2_tolerance=ms2_tolerance,
        calibration_set_size=calibration_set_size,
        rescoring_engine=rescoring_engine,
        rng=rng,
        test_fdr=test_fdr,
        spectrum_path=spectrum_path,
        fasta_file=fasta_file,
        id_decoy_pattern=id_decoy_pattern,
        lower_score_is_better=lower_score_is_better,
        log_level=log_level,
        spectrum_id_pattern=spectrum_id_pattern,
        psm_id_pattern=psm_id_pattern
    )
    logging.info("MS²Rescore config:")
    logging.info(config)
    rescore_idxml(psm_file, output_path, config)


def main(**kwargs):
    config = parse_cli_arguments_to_config(**kwargs)
    rescore_idxml(kwargs["psm_file"], kwargs["output_path"], config)
