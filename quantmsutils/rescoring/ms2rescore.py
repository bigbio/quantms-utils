# Written by Jonas Scheid under the MIT license
# Contributions by Yasset Perez-Riverol and Dai Chengxin
# This script is part of the quantmsutils package

import importlib.resources
import json
import logging

import click
import pyopenms as oms
from ms2rescore import package_data, rescore
from psm_utils import PSMList
from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
from typing import Iterable, List, Union
from pathlib import Path
from psm_utils.psm import PSM

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


class IDXMLReaderPatch(IdXMLReader):
    def __init__(self, filename: Union[Path, str], *args, **kwargs) -> None:
        """
        Patch Reader for idXML files based on IDXMLReader.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to idXML file.

        Examples
        --------
        """
        super().__init__(filename, *args, **kwargs)
        self.protein_ids, self.peptide_ids = self._parse_idxml()
        self.user_params_metadata = self._get_userparams_metadata(self.peptide_ids[0].getHits()[0])
        self.rescoring_features = self._get_rescoring_features(self.peptide_ids[0].getHits()[0])
        self.skip_invalid_psm = 0

    def __iter__(self) -> Iterable[PSM]:
        """
        Iterate over file and return PSMs one-by-one.
                Test cases will:

        Input PSM 1: PeptideHit	with metavalue
            "MSGF:ScoreRatio" value="0.212121212121212"/>
            "MSGF:Energy" value="130.0"/>
            "MSGF:lnEValue" value="-3.603969939390662"/>
            "MSGF:lnExplainedIonCurrentRatio" value="-0.881402756873971"/>
            "MSGF:lnNTermIonCurrentRatio" value="-1.931878317286471"/>
            "MSGF:lnCTermIonCurrentRatio" value="-1.311462733724937"/>
            "MSGF:lnMS2IonCurrent" value="9.702930189540499"/>
            "MSGF:MeanErrorTop7" value="259.986879999999985"/>
            "MSGF:sqMeanErrorTop7" value="6.75931777721344e04"/>
            "MSGF:StdevErrorTop7" value="143.678020000000004"/>
        PSM2: PeptideHit No above metaValue

        Run:
        reader = IDXMLReaderPatch(input_file)
        psm_list = reader.read_file()

        psm_list: return [PSM 1]

        """
        for peptide_id in self.peptide_ids:
            for peptide_hit in peptide_id.getHits():
                psm = self._parse_psm(self.protein_ids, peptide_id, peptide_hit)
                if psm is not None:
                    yield psm
                else:
                    self.skip_invalid_psm += 1

    def _parse_psm(
            self,
            protein_ids: oms.ProteinIdentification,
            peptide_id: oms.PeptideIdentification,
            peptide_hit: oms.PeptideHit,
    ) -> PSM:
        """
        Parse idXML :py:class:`~pyopenms.PeptideHit` to :py:class:`~psm_utils.psm.PSM`.

        Uses additional information from :py:class:`~pyopenms.ProteinIdentification` and
        :py:class:`~pyopenms.PeptideIdentification` to annotate parameters of the
        :py:class:`~psm_utils.psm.PSM` object.
        """
        peptidoform = self._parse_peptidoform(
            peptide_hit.getSequence().toString(), peptide_hit.getCharge()
        )
        # This is needed to calculate a qvalue before rescoring the PSMList
        peptide_id_metadata = {
            "idxml:score_type": str(peptide_id.getScoreType()),
            "idxml:higher_score_better": str(peptide_id.isHigherScoreBetter()),
            "idxml:significance_threshold": str(peptide_id.getSignificanceThreshold()),
        }
        peptide_hit_metadata = {
            key: peptide_hit.getMetaValue(key) for key in self.user_params_metadata
        }

        # Get search engines score features and check valueExits
        rescoring_features = {}
        for key in self.rescoring_features:
            feature = peptide_hit.metaValueExists(key)
            if not feature:
                return None
            else:
                rescoring_features[key] = float(peptide_hit.getMetaValue(key))

        return PSM(
            peptidoform=peptidoform,
            spectrum_id=peptide_id.getMetaValue("spectrum_reference"),
            run=self._get_run(protein_ids, peptide_id),
            is_decoy=self._is_decoy(peptide_hit),
            score=peptide_hit.getScore(),
            precursor_mz=peptide_id.getMZ(),
            retention_time=peptide_id.getRT(),
            # NOTE: ion mobility will be supported by OpenMS in the future
            protein_list=[
                accession.decode() for accession in peptide_hit.extractProteinAccessionsSet()
            ],
            rank=peptide_hit.getRank() + 1,  # 0-based to 1-based
            source="idXML",
            # Storing proforma notation of peptidoform and UNIMOD peptide sequence for mapping back
            # to original sequence in writer
            provenance_data={str(peptidoform): peptide_hit.getSequence().toString()},
            # Store metadata of PeptideIdentification and PeptideHit objects
            metadata={**peptide_id_metadata, **peptide_hit_metadata},

            rescoring_features=rescoring_features,
        )


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
    reader = IDXMLReaderPatch(input_file)
    psm_list = reader.read_file()

    if reader.skip_invalid_psm != 0:
        logging.warning(
            f"Removed {reader.skip_invalid_psm} PSMs without search engine features!"
        )

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
