import pytest
from click.testing import CliRunner

from pyquantms.pyquantmsc import cli


def test_future_test():
    assert 1 == 1

# test for the create_diann_cfg command in cli
def test_create_diann_cfg():
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['create_diann_cfg','--help'])

    assert result.exit_code == 0