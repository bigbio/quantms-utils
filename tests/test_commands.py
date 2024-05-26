from click.testing import CliRunner
from pyquantms.pyquantmsc import cli


def test_future_test():
    assert 1 == 1


# test for the create_diann_cfg command in cli
def test_create_diann_cfg():
    runner = CliRunner()
    result = runner.invoke(cli, ["dianncfg", "--help"])

    assert result.exit_code == 0


# test for the add_sage_feature command in cli
def test_add_sage_feature():
    runner = CliRunner()
    result = runner.invoke(cli, ["add_sage_feature", "--help"])

    assert result.exit_code == 0


# test for the mzml_statistics command in cli
def test_mzml_statistics():
    runner = CliRunner()
    result = runner.invoke(cli, ["mzml_statistics", "--help"])

    assert result.exit_code == 0


# test for the diann_convert command in cli
def test_diann_convert():
    runner = CliRunner()
    result = runner.invoke(cli, ["diann_convert", "--help"])

    assert result.exit_code == 0


# test for the extract_sample_from_expdesign command in cli
def test_extract_sample_from_expdesign():
    runner = CliRunner()
    result = runner.invoke(cli, ["extract_sample", "--help"])

    assert result.exit_code == 0


# test for the rescoring command in cli
def test_ms2rescore():
    runner = CliRunner()
    result = runner.invoke(cli, ["ms2rescore", "--help"])

    assert result.exit_code == 0


# test for the convert_psm command in cli
def test_convert_psm():
    runner = CliRunner()
    result = runner.invoke(cli, ["convert_psm", "--help"])

    assert result.exit_code == 0


# test for the check_samplesheet command in cli
def test_check_samplesheet():
    runner = CliRunner()
    result = runner.invoke(cli, ["check_samplesheet", "--help"])

    assert result.exit_code == 0

# test the validation of an SDRF file
def test_check_samplesheet_sdrf():
    runner = CliRunner()
    result = runner.invoke(cli, ["check_samplesheet", "--is_sdrf", "--check_ms", "--input", "tests/test_data/PXD000001.sdrf.tsv"])

    assert result.exit_code == 0

# test extract_sample_from_expdesign command in cli
def test_extract_sample_from_expdesign():
    runner = CliRunner()
    result = runner.invoke(cli, ["extract_sample", "--expdesign", "tests/test_data/BSA_design_urls.tsv"])

    assert result.exit_code == 0

# test psm conversion command in cli
def test_convert_psm():
    runner = CliRunner()
    result = runner.invoke(cli, ["convert_psm", "--idxml", "tests/test_data/BSA1_F1_msgf_idx_fdr_idpep_switched_filter.idXML", "--spectra_file", "tests/test_data/BSA1_F1_spectrum_df.csv"])

    assert result.exit_code == 0

# test mzml statistics command in cli
def test_mzml_statistics():
    runner = CliRunner()
    result = runner.invoke(cli, ["mzml_statistics", "--ms_path", "tests/test_data/BSA1_F1.mzML"])

    assert result.exit_code == 0