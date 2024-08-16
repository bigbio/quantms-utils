import filecmp

from click.testing import CliRunner

from quantmsutils.quantmsutilsc import cli


# test for the create_diann_cfg command in cli
def test_create_diann_cfg_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["dianncfg", "--help"])

    assert result.exit_code == 0


# test for the add_sage_feature command in cli
def test_add_sage_feature_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["sage2feature", "--help"])

    assert result.exit_code == 0


# test for the mzml_statistics command in cli
def test_mzml_statistics_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["mzmlstats", "--help"])

    assert result.exit_code == 0


# test for the diann_convert command in cli
def test_diann_convert_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["diann2mztab", "--help"])

    assert result.exit_code == 0


# test for the extract_sample_from_expdesign command in cli
def test_extract_sample_from_expdesign_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["openms2sample", "--help"])

    assert result.exit_code == 0


# test for the rescoring command in cli
def test_ms2rescore_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["ms2rescore", "--help"])

    assert result.exit_code == 0


def test_ms2rescore():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "ms2rescore",
            "--psm_file",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML",
            "--spectrum_path",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML",
            "--processes",
            "2",
            "--ms2pip_model",
            "HCD2021",
            "--feature_generators",
            "'ms2pip,deeplc'",
            "--id_decoy_pattern",
            "^rev",
            "--test_fdr",
            "0.05",
        ],
    )
    assert result.exit_code == 0


def test_sage_feature_file():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "sage2feature",
            "--idx_file",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore.idXML",
            "--output_file",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore_feat_gen.idXML",
            "--feat_file",
            "tests/test_data/tmt_erwinia_1ulsike_top10hcd_isol2_45stepped_60min_01_sage_ms2rescore.idxml.feature_names.tsv",
        ],
    )

    assert result.exit_code == 0


def test_spectrum2fature_file():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "spectrum2feature",
            "--ms_path",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML",
            "--idxml",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore.idXML",
            "--output",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore_snr.idXML",
        ],
    )

    assert result.exit_code == 0


def test_diann2mztab_example():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "diann2mztab",
            "--folder",
            "tests/test_data/diann2mztab/",
            "--exp_design",
            "tests/test_data/diann2mztab/PXD026600.sdrf_openms_design.tsv",
            "--diann_version",
            "tests/test_data/diann2mztab/versions.yml",
            "--dia_params",
            "20.0;ppm;10.0;ppm;Trypsin;Carbamidomethyl (C);Oxidation (M)",
            "--charge",
            "3",
            "--missed_cleavages",
            "1",
            "--qvalue_threshold",
            "0.01",
        ],
    )
    assert result.exit_code == 0


# test for the convert_psm command in cli
def test_convert_psm_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["psmconvert", "--help"])

    assert result.exit_code == 0


# test for the check_samplesheet command in cli
def test_check_samplesheet_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["checksamplesheet", "--help"])

    assert result.exit_code == 0


# test the validation of an SDRF file
def test_check_samplesheet_sdrf():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "checksamplesheet",
            "--is_sdrf",
            "--exp_design",
            "tests/test_data/PXD000001.sdrf.tsv",
        ],
    )

    assert result.exit_code == 0


# test extract_sample_from_expdesign command in cli
def test_extract_sample_from_expdesign():
    runner = CliRunner()
    result = runner.invoke(
        cli, ["openms2sample", "--expdesign", "tests/test_data/BSA_design_urls.tsv"]
    )

    assert result.exit_code == 0


# test psm conversion command in cli
def test_convert_psm():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "psmconvert",
            "--idxml",
            "tests/test_data/BSA1_F1_msgf_idx_fdr_idpep_switched_filter.idXML",
            "--spectra_file",
            "tests/test_data/BSA1_F1_spectrum_df.parquet",
        ],
    )

    assert result.exit_code == 0


# test mzml statistics command in cli
def test_mzml_statistics():
    runner = CliRunner()
    result = runner.invoke(
        cli, ["mzmlstats", "--ms_path", "tests/test_data/BSA1_F1.mzML"]
    )

    assert result.exit_code == 0
