from pathlib import Path
import os
import pandas as pd
import pytest
from click.testing import CliRunner

from quantmsutils.mzml.ms1_feature_finder import MS1FeatureDetector
from quantmsutils.quantmsutilsc import cli

# Define constants at the top for better maintainability
TESTS_DIR = Path(__file__).parent
TEST_DATA_DIR = TESTS_DIR / "test_data"
DIANN_TEST_DIR = TEST_DATA_DIR / "diann2mztab"

TMT_MZML_FILE = TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
TMT_MS1_FEAURES = (
    TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_ms1_features.parquet"
)
TMT_IDXML_FILE = (
    TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML"
)
TMT_MS_INFO_FILE = (
    TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_ms_info.parquet"
)
TMT_STATIC_MS2_INFO_FILE = (
    TEST_DATA_DIR / "STATIC_TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_ms2_info.parquet"
)
TMT_MS2_FILE = (
    TEST_DATA_DIR / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_ms2_info.parquet"
)


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
            "--folder",
            str(DIANN_TEST_DIR),
            "--exp_design",
            str(DIANN_TEST_DIR / "PXD026600.sdrf_openms_design.tsv"),
            "--diann_version",
            str(DIANN_TEST_DIR / "versions.yml"),
            "--dia_params",
            "20.0;ppm;10.0;ppm;Trypsin;Carbamidomethyl (C);Oxidation (M)",
            "--charge",
            "3",
            "--missed_cleavages",
            "1",
            "--qvalue_threshold",
            "0.01",
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
            "--exp_design",
            str(TEST_DATA_DIR / "PXD000001.sdrf.tsv"),
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
    """
    Test class for PSM conversion commands, it takes an ms2 file in parquet and an idXML file
    and converts it to a parquet file with PSM information.
    """

    def test_convert_psm(self):
        """Test converting PSM data"""
        args = [
            "--idxml",
            str(TMT_IDXML_FILE),
            "--ms2_file",
            str(TMT_STATIC_MS2_INFO_FILE),
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
        if TMT_MS_INFO_FILE.exists():
            TMT_MS_INFO_FILE.unlink()
        if TMT_MS2_FILE.exists():
            TMT_MS2_FILE.unlink()
        yield
        # Teardown - can add cleanup code here if needed

    def test_mzml_statistics(self):
        """Test mzML statistics on BSA sample"""
        args = ["--ms2_file", "--ms_path", str(TMT_MZML_FILE)]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0
        assert TMT_MS_INFO_FILE.exists(), "Output file was not created"

        # Compare with reference data
        output_table = pd.read_parquet(TMT_MS_INFO_FILE)

        assert len(output_table) > 0

        expected_columns = ["scan", "rt", "ms_level", "precursor_charge"]
        for col in expected_columns:
            assert col in output_table.columns, f"Expected column {col} missing from output"

    @pytest.mark.skip("Test to be run locally, with big files")
    def test_mzml_statistics_local(self):

        args = [
            "--ms2_file",
            "--ms_path",
            str(TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1.mzML"),
        ]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0
        assert os.path.exists(TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1_ms_info.parquet")

        output_table = pd.read_parquet(
            TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1_ms_info.parquet"
        )
        assert len(output_table) > 0, "Output table is empty"

    @pytest.mark.skip("Test to be run locally, with bruker file")
    def test_mzml_statistics_bruker(self):
        """Test mzML statistics on Bruker sample"""
        args = [
            "--ms2_file",
            "--ms_path",
            str(TEST_DATA_DIR / "hMICAL1_coiPAnP-N2-200_3Murea-1Mthiourea-200mMtcep_14733.d"),
        ]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0
        assert os.path.exists(
            TEST_DATA_DIR
            / "hMICAL1_coiPAnP-N2-200_3Murea-1Mthiourea-200mMtcep_14733_ms_info.parquet"
        )

        output_table = pd.read_parquet(
            TEST_DATA_DIR
            / "hMICAL1_coiPAnP-N2-200_3Murea-1Mthiourea-200mMtcep_14733_ms_info.parquet"
        )
        assert len(output_table) > 0, "Output table is empty"


class TestFeatureFinder:

    def test_feature_finder(self):
        """Test feature finder on TMT data"""

        detector = MS1FeatureDetector(min_ptic=0.05, max_ptic=0.95)
        result = detector.process_file(input_file=TMT_MZML_FILE, output_file=TMT_MS1_FEAURES)

        parquet_df = pd.read_parquet(result)
        assert not parquet_df.empty, "Output table is empty"

        # result_dia = detector.process_file(
        #     input_file=TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1.mzML",
        #     output_file=TEST_DATA_DIR / "RD139_Narrow_UPS1_0_1fmol_inj1_ms1_features.parquet",
        # )
        #
        # parquet_df_dia = pd.read_parquet(result_dia)
        # assert not parquet_df_dia.empty, "Output table is empty"
