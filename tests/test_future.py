from click.testing import CliRunner

from pyquantms.pyquantmsc import cli


def test_future_test():
    assert 1 == 1

# test for the create_diann_cfg command in cli
def test_create_diann_cfg():
    runner = CliRunner()
    result = runner.invoke(cli,['create_diann_cfg','--help'])

    assert result.exit_code == 0

# test for the add_sage_feature command in cli
def test_add_sage_feature():
    runner = CliRunner()
    result = runner.invoke(cli,['add_sage_feature','--help'])

    assert result.exit_code == 0

# test for the mzml_statistics command in cli
def test_mzml_statistics():
    runner = CliRunner()
    result = runner.invoke(cli,['mzml_statistics','--help'])

    assert result.exit_code == 0

# test for the diann_convert command in cli
def test_diann_convert():
    runner = CliRunner()
    result = runner.invoke(cli,['diann_convert','--help'])

    assert result.exit_code == 0

# test for the extract_sample_from_expdesign command in cli
def test_extract_sample_from_expdesign():
    runner = CliRunner()
    result = runner.invoke(cli,['extract_sample','--help'])

    assert result.exit_code == 0