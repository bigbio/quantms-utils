from pathlib import Path
import os
import pandas as pd
import pytest
from click.testing import CliRunner
from quantmsutils.quantmsutilsc import cli

# Define constants at the top for better maintainability
TESTS_DIR = Path(__file__).parent
TEST_DATA_DIR = TESTS_DIR / "test_data"
DIANN_TEST_DIR = TEST_DATA_DIR / "diann2mztab"
BSA_IDXML_FILE = TEST_DATA_DIR / "BSA1_F1_msgf_idx_fdr_idpep_switched_filter.idXML"
BSA_SPECTRA_FILE = TEST_DATA_DIR / "BSA1_F1_spectrum_df.parquet"
BSA_MZML_FILE = TEST_DATA_DIR / "BSA1_F1.mzML"
BSA_MS_INFO_FILE = TEST_DATA_DIR / "BSA1_F1_ms_info.parquet"
BSA_TEST_MS_INFO_FILE = TEST_DATA_DIR / "BSA1_F1_test_ms_info.parquet"
TMT_MZML_FILE = TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
TMT_MS_INFO_FILE = TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_ms_info.parquet"


# Helper function to create a CLI runner and run a command
def run_cli_command(command, args=None):
    runner = CliRunner()
    if args:
        result = runner.invoke(cli, [command] + args)
    else:
        result = runner.invoke(cli, [command, "--help"])
    return result


class TestCLIHelpMessages:
    """Test class for CLI help messages"""

    @pytest.mark.parametrize(
        "command",
        [
            "dianncfg",
            "mzmlstats",
            "diann2mztab",
            "openms2sample",
            "psmconvert",
            "checksamplesheet",
        ],
    )
    def test_help_messages(self, command):
        """Test all CLI help messages with a parametrized test"""
        result = run_cli_command(command)
        assert result.exit_code == 0
        assert "Usage:" in result.output


class TestDiannCommands:
    """Test class for DIA-NN related commands"""

    def test_diann2mztab_example(self):
        """Test the DIA-NN to mzTab conversion with example data"""
        args = [
            "--folder", str(DIANN_TEST_DIR),
            "--exp_design", str(DIANN_TEST_DIR / "PXD026600.sdrf_openms_design.tsv"),
            "--diann_version", str(DIANN_TEST_DIR / "versions.yml"),
            "--dia_params", "20.0;ppm;10.0;ppm;Trypsin;Carbamidomethyl (C);Oxidation (M)",
            "--charge", "3",
            "--missed_cleavages", "1",
            "--qvalue_threshold", "0.01",
        ]
        result = run_cli_command("diann2mztab", args)
        assert result.exit_code == 0
        # Additional assertions could check for expected output files


class TestSamplesheetCommands:
    """Test class for samplesheet related commands"""

    def test_check_samplesheet_sdrf(self):
        """Test the validation of an SDRF file"""
        args = [
            "--is_sdrf",
            "--exp_design", str(TEST_DATA_DIR / "PXD000001.sdrf.tsv"),
        ]
        result = run_cli_command("checksamplesheet", args)
        assert result.exit_code == 0

    def test_extract_sample_from_expdesign(self):
        """Test extracting sample information from experiment design"""
        args = ["--expdesign", str(TEST_DATA_DIR / "BSA_design_urls.tsv")]
        result = run_cli_command("openms2sample", args)
        assert result.exit_code == 0
        # Could add assertions to check for expected output files


class TestPSMConversion:
    """Test class for PSM conversion commands"""

    def test_convert_psm(self):
        """Test converting PSM data"""
        args = [
            "--idxml", str(BSA_IDXML_FILE),
            "--spectra_file", str(BSA_SPECTRA_FILE),
        ]
        result = run_cli_command("psmconvert", args)
        assert result.exit_code == 0
        # Could add assertions to check for expected output files


class TestMzMLStatistics:
    """Test class for mzML statistics commands"""

    @pytest.fixture(autouse=True)
    def setup_teardown(self):
        """Setup and teardown for mzML statistics tests"""
        # Setup - remove output files if they exist
        if BSA_MS_INFO_FILE.exists():
            BSA_MS_INFO_FILE.unlink()
        if TMT_MS_INFO_FILE.exists():
            TMT_MS_INFO_FILE.unlink()
        yield
        # Teardown - can add cleanup code here if needed

    def test_mzml_statistics_bsa(self):
        """Test mzML statistics on BSA sample"""
        args = ["--id_only", "--ms_path", str(BSA_MZML_FILE)]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0
        assert BSA_MS_INFO_FILE.exists(), "Output file was not created"

        # Compare with reference data
        output_table = pd.read_parquet(BSA_MS_INFO_FILE).set_index("scan")
        reference_table = pd.read_parquet(BSA_TEST_MS_INFO_FILE).set_index("scan")

        assert len(output_table) == len(reference_table), "Number of scans doesn't match"
        # Could add more specific comparisons of table contents

    def test_mzml_statistics_tmt(self):
        """Test mzML statistics on TMT sample"""
        args = ["--id_only", "--ms_path", str(TMT_MZML_FILE)]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0
        assert TMT_MS_INFO_FILE.exists(), "Output file was not created"

        output_table = pd.read_parquet(TMT_MS_INFO_FILE)
        assert len(output_table) > 0, "Output table is empty"

        # Could add assertions for expected columns or data validation
        expected_columns = ["scan", "rt", "ms_level", "precursor_charge"]
        for col in expected_columns:
            assert col in output_table.columns, f"Expected column {col} missing from output"

    @pytest.mark.skip("Test to be run locally, with big files")
    def test_mzml_statistics_local(self):

        args = ["--id_only", "--ms_path", str(TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1.mzML")]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0
        assert os.path.exists(TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1_ms_info.parquet")

        output_table = pd.read_parquet(TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1_ms_info.parquet")
        assert len(output_table) > 0, "Output table is empty"