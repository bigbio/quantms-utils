"""Regression tests for lossless, design-aware LFQ partitioning."""

import pytest
from click.testing import CliRunner

from quantmsutils.quantmsutilsc import cli
from quantmsutils.sdrf import split_lfq_groups as groups

COLUMNS = ["comment[instrument]", "comment[gradient duration]"]


@pytest.fixture
def sdrf():
    return [
        [
            "source name",
            "comment[data file]",
            *COLUMNS,
            "comment[modification parameters]",
            "comment[modification parameters]",
        ],
        ["sample A", "a.raw", "instrument A", "30 min", "fixed C", "variable M"],
        ["sample B", "b.raw", "instrument B", "120 min", "fixed C", "variable M"],
        ["sample C", "c.raw", "instrument A", "30 min", "fixed C", "variable M"],
        ["sample D", "d.raw", "instrument B", "120 min", "fixed C", "variable M"],
    ]


@pytest.fixture
def design():
    return [
        ["Fraction_Group", "Fraction", "Spectra_Filepath", "Label", "Sample"],
        ["1", "1", "a.mzML", "1", "1"],
        ["2", "1", "b.mzML", "1", "2"],
        ["3", "1", "c.mzML", "1", "3"],
        ["4", "1", "d.mzML", "1", "4"],
        [],
        ["Sample", "MSstats_Condition", "MSstats_BioReplicate"],
        ["1", "control", "D1"],
        ["2", "control", "D2"],
        ["3", "spike", "D3"],
        ["4", "spike", "D4"],
    ]


def write_inputs(root, sdrf, design):
    sdrf_path, design_path = root / "full.sdrf.tsv", root / "full_design.tsv"
    groups.write_table(sdrf_path, sdrf)
    groups.write_table(design_path, design)
    return sdrf_path, design_path


def run_split(root, sdrf, design, columns=COLUMNS):
    sdrf_path, design_path = write_inputs(root, sdrf, design)
    before = sdrf_path.read_bytes(), design_path.read_bytes()
    result = groups.split_groups(sdrf_path, design_path, columns, root / "groups")
    assert before == (sdrf_path.read_bytes(), design_path.read_bytes())
    return result


def test_complete_partition_and_replicate_identity(tmp_path, sdrf, design):
    manifest = run_split(tmp_path, sdrf, design)
    assert {row[1] for row in manifest[1:]} == {"a", "b", "c", "d"}
    assert len(manifest) == 5
    assert len({row[0] for row in manifest[1:]}) == 2
    emitted_rows = []
    seen_replicates = set()
    for name in {row[2] for row in manifest[1:]}:
        rows = groups.read_table(tmp_path / "groups" / name)
        assert rows[0] == sdrf[0]
        emitted_rows.extend(rows[1:])
    for name in {row[3] for row in manifest[1:]}:
        files, samples = groups.design_tables(groups.read_table(tmp_path / "groups" / name))
        assert {row[4] for row in files[1:]} == {row[0] for row in samples[1:]}
        assert {row[1] for row in samples[1:]} == {"control", "spike"}
        seen_replicates.update(row[2] for row in samples[1:])
    assert sorted(emitted_rows) == sorted(sdrf[1:])
    assert seen_replicates == {"D1", "D2", "D3", "D4"}


def test_group_ids_survive_row_reordering(tmp_path, sdrf, design):
    first = run_split(tmp_path, sdrf, design)
    sdrf[1:] = reversed(sdrf[1:])
    second = run_split(tmp_path, sdrf, design)
    assert {tuple(row[:2]) for row in first[1:]} == {tuple(row[:2]) for row in second[1:]}


def test_unknown_grouping_column(tmp_path, sdrf, design):
    with pytest.raises(ValueError, match="Expected exactly one"):
        run_split(tmp_path, sdrf, design, columns=["missing"])


def test_missing_group_value(tmp_path, sdrf, design):
    sdrf[1][2] = "not available"
    with pytest.raises(ValueError, match="Missing LFQ grouping"):
        run_split(tmp_path, sdrf, design)


def test_conflicting_run_assignment(tmp_path, sdrf, design):
    sdrf.append(["sample A", "a.raw", "instrument B", "120 min", "fixed C", "variable M"])
    with pytest.raises(ValueError, match="multiple LFQ groups"):
        run_split(tmp_path, sdrf, design)


def test_unmatched_run_is_not_dropped(tmp_path, sdrf, design):
    sdrf[1][1] = "unexpected.raw"
    with pytest.raises(ValueError, match="run sets differ"):
        run_split(tmp_path, sdrf, design)


def test_fractionated_sample_cannot_cross_groups(tmp_path, sdrf, design):
    design[2][0] = "1"
    with pytest.raises(ValueError, match="crosses LFQ groups"):
        run_split(tmp_path, sdrf, design)


def test_missing_sample_metadata_is_rejected(tmp_path, sdrf, design):
    design[-1][0] = "5"
    with pytest.raises(ValueError, match="missing sample metadata"):
        run_split(tmp_path, sdrf, design)


def test_cli_takes_one_column_option_per_grouping_column(tmp_path, sdrf, design):
    sdrf_path, design_path = write_inputs(tmp_path, sdrf, design)
    output = tmp_path / "groups"
    arguments = ["splitlfqgroups", "--sdrf", str(sdrf_path), "--design", str(design_path)]
    for column in COLUMNS:
        arguments += ["--column", column]
    result = CliRunner().invoke(cli, [*arguments, "--output", str(output)])

    assert result.exit_code == 0, result.output
    manifest = groups.read_table(output / "lfq_groups.tsv")
    assert manifest[0] == ["group_id", "run_id", "sdrf_file", "design_file", *COLUMNS]
    assert sorted(row[1] for row in manifest[1:]) == ["a", "b", "c", "d"]
