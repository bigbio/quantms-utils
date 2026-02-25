from pathlib import Path
import os
import tempfile
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from quantmsutils.mzml.ms1_feature_finder import MS1FeatureDetector
from quantmsutils.psm.psm_conversion import mods_position
from quantmsutils.quantmsutilsc import cli

# Define constants at the top for better maintainability
TESTS_DIR = Path(__file__).parent
TEST_DATA_DIR = TESTS_DIR / "test_data"
DIANN_TEST_DIR = TEST_DATA_DIR / "diann2msstats"

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
BRUKER_D_FILE = TEST_DATA_DIR / "hMICAL1_coiPAnP-N2-200_3Murea-1Mthiourea-200mMtcep_14733.d"
BRUKER_MS_INFO_FILE = (
    TEST_DATA_DIR / "hMICAL1_coiPAnP-N2-200_3Murea-1Mthiourea-200mMtcep_14733_ms_info.parquet"
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
            "diann2msstats",
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

    def test_diann2msstats_example(self):
        """Test the DIA-NN to MSstats conversion with example data"""
        report_path = (DIANN_TEST_DIR / "diann_report.tsv").resolve()
        exp_design_path = (DIANN_TEST_DIR / "PXD026600.sdrf_openms_design.tsv").resolve()
        assert report_path.exists(), f"Test report missing: {report_path}"
        assert exp_design_path.exists(), f"Test design missing: {exp_design_path}"

        args = [
            "--report",
            str(report_path),
            "--exp_design",
            str(exp_design_path),
            "--qvalue_threshold",
            "0.01",
        ]
        result = run_cli_command("diann2msstats", args)
        if result.exit_code != 0:
            raise AssertionError(
                f"diann2msstats failed (exit {result.exit_code}). "
                f"stdout: {result.output!r}, stderr: {result.stderr!r}"
            )

    def test_dianncfg_example(self):
        """Test generating the DIA-NN config with example data"""
        args = [
            "--enzyme",
            "Trypsin",
            "--fix_mod",
            "Carbamidomethyl (C)",
            "--var_mod",
            "Oxidation (M),Phospho (S),Phospho (T),Phospho (Y),Acetyl (Protein N-term),Acetyl (K),Acetyl (R),Met-loss (Protein N-term M)",
        ]
        result = run_cli_command("dianncfg", args)

        assert result.exit_code == 0


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
        """Test converting PSM data with ms2_file"""
        args = [
            "--idxml",
            str(TMT_IDXML_FILE),
            "--ms2_file",
            str(TMT_STATIC_MS2_INFO_FILE),
        ]
        result = run_cli_command("psmconvert", args)
        assert result.exit_code == 0

        output_file = Path(TMT_IDXML_FILE).with_name(
            f"{Path(TMT_IDXML_FILE).stem}_psm.parquet"
        )
        assert output_file.exists(), f"PSM output file was not created: {output_file}"
        df = pd.read_parquet(output_file)
        assert len(df) > 0, "PSM output parquet is empty"
        assert "sequence" in df.columns
        assert "scan_number" in df.columns

    def test_convert_psm_without_ms2(self):
        """Test converting PSM data without ms2_file (regression test for None dereference)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "test_psm.parquet")
            args = [
                "--idxml",
                str(TMT_IDXML_FILE),
                "--output_file",
                output_file,
            ]
            result = run_cli_command("psmconvert", args)
            assert result.exit_code == 0

            df = pd.read_parquet(output_file)
            assert len(df) > 0, "PSM output parquet is empty"
            assert "sequence" in df.columns


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

    @pytest.mark.skipif(
        not BRUKER_D_FILE.exists(),
        reason="Bruker .d test data not available (not tracked in git)",
    )
    def test_mzml_statistics_bruker(self):
        """Test mzML statistics on Bruker .d file (regression test for if/elif bug)"""
        if BRUKER_MS_INFO_FILE.exists():
            BRUKER_MS_INFO_FILE.unlink()

        args = [
            "--ms_path",
            str(BRUKER_D_FILE),
        ]
        result = run_cli_command("mzmlstats", args)

        assert result.exit_code == 0, f"Bruker .d processing failed: {result.output}"
        assert BRUKER_MS_INFO_FILE.exists(), "Bruker output parquet was not created"

        output_table = pd.read_parquet(BRUKER_MS_INFO_FILE)
        assert len(output_table) > 0, "Output table is empty"

        expected_columns = ["scan", "ms_level"]
        for col in expected_columns:
            assert col in output_table.columns, f"Expected column {col} missing from Bruker output"


class TestFeatureFinder:

    def test_feature_finder(self):
        """Test feature finder on TMT data"""

        detector = MS1FeatureDetector()
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


class TestModsPosition:
    """Unit tests for the mods_position function in psm_conversion"""

    def test_single_modification(self):
        """Test peptide with a single modification"""
        result = mods_position("PEPTM(Oxidation)IDE")
        assert result == ["5-Oxidation"]

    def test_multiple_modifications(self):
        """Test peptide with multiple modifications"""
        result = mods_position("PEC(Carbamidomethyl)PTMC(Carbamidomethyl)IDE")
        assert result == ["3-Carbamidomethyl", "7-Carbamidomethyl"]

    def test_no_modifications(self):
        """Test unmodified peptide returns NaN"""
        result = mods_position("PEPTIDE")
        assert result is np.nan or (isinstance(result, float) and np.isnan(result))

    def test_leading_dot(self):
        """Test peptide with leading dot (sometimes from OpenMS output)"""
        result = mods_position(".PEPTM(Oxidation)IDE")
        assert result == ["5-Oxidation"]

    def test_nterm_modification(self):
        """Test N-terminal modification"""
        result = mods_position("(Acetyl)PEPTIDE")
        assert result == ["0-Acetyl"]


class TestExtractSampleMixture:
    """Test extract_sample with MSstats_Mixture column (covers DataFrame.append fix)"""

    def test_extract_sample_with_mixture(self):
        """Test that extract_sample works with MSstats_Mixture in the experiment design"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a test experiment design file with MSstats_Mixture
            design_file = os.path.join(tmpdir, "test_design.tsv")
            with open(design_file, "w") as f:
                f.write("Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n")
                f.write("1\t1\tfile1.mzML\t1\t1\n")
                f.write("2\t1\tfile2.mzML\t1\t2\n")
                f.write("\n")
                f.write("Sample\tMSstats_Condition\tMSstats_BioReplicate\tMSstats_Mixture\n")
                f.write("1\tCondition_A\t1\tMixture_1\n")
                f.write("2\tCondition_B\t2\tMixture_1\n")

            # Run extract_sample from within the temp dir so output goes there
            original_dir = os.getcwd()
            try:
                os.chdir(tmpdir)
                args = ["--expdesign", design_file]
                result = run_cli_command("openms2sample", args)
                assert result.exit_code == 0

                output_file = os.path.join(tmpdir, "test_design_sample.csv")
                assert os.path.exists(output_file), "Sample output file was not created"

                df = pd.read_csv(output_file, sep="\t")
                assert len(df) == 2, f"Expected 2 rows, got {len(df)}"
                assert "Spectra_Filepath" in df.columns
                assert "Sample" in df.columns
                # The Sample column should contain mixture IDs
                assert df["Sample"].iloc[0] == "Mixture_1"
            finally:
                os.chdir(original_dir)
