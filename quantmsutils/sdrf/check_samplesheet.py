import logging
import sys

import click
import pandas as pd

from sdrf_pipelines.sdrf.sdrf import read_sdrf

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Minimal columns required to run quantms/quantmsdiann pipelines.
# These are checked in --minimal mode instead of full schema validation.
MINIMAL_REQUIRED_COLUMNS = [
    "source name",
    "assay name",
    "comment[data file]",
    "comment[label]",
    "comment[cleavage agent details]",
    "comment[instrument]",
    "comment[proteomics data acquisition method]",
    "technology type",
]

# Recommended columns: warn if missing but don't fail
MINIMAL_RECOMMENDED_COLUMNS = [
    "comment[precursor mass tolerance]",
    "comment[fragment mass tolerance]",
    "comment[dissociation method]",
    "comment[technical replicate]",
    "comment[fraction identifier]",
]


def check_sdrf(
    input_sdrf: str,
    template: str = "ms-proteomics",
    minimal: bool = False,
    use_ols_cache_only: bool = False,
):
    """
    Check the SDRF file for errors.

    :param input_sdrf: Path to the SDRF file to check
    :param template: Schema template for full validation (e.g. 'ms-proteomics', 'dia-acquisition')
    :param minimal: Only validate columns required to run the pipeline (skip organism, etc.)
    :param use_ols_cache_only: Use OLS cache instead of live OLS service
    """
    if minimal:
        errors = _validate_minimal(input_sdrf)
    else:
        df = read_sdrf(input_sdrf)
        errors = df.validate_sdrf(
            template=template,
            use_ols_cache_only=use_ols_cache_only,
        )

    for error in errors:
        print(error)

    sys.exit(bool(errors))


def _validate_minimal(input_sdrf: str) -> list[str]:
    """Validate only the columns required to run the pipeline.

    Returns a list of error strings. Only missing required columns
    produce errors; missing recommended columns produce warnings (non-blocking).
    """
    df_header = pd.read_csv(input_sdrf, sep="\t", nrows=0)
    columns_lower = [c.lower() for c in df_header.columns]
    errors = []

    # Reject header-only files
    df_rows = pd.read_csv(input_sdrf, sep="\t", nrows=1)
    if len(df_rows) == 0:
        errors.append("ERROR: SDRF file contains a header but no data rows.")
        return errors

    # Check required columns (case-insensitive)
    for col in MINIMAL_REQUIRED_COLUMNS:
        if col.lower() not in columns_lower:
            errors.append(f"ERROR: Required column '{col}' is missing from the SDRF file.")

    # Check at least one modification parameters column exists
    has_mod_col = any(c.startswith("comment[modification parameters") for c in columns_lower)
    if not has_mod_col:
        errors.append(
            "ERROR: At least one 'comment[modification parameters]' column is required."
        )

    # Warn about recommended columns (non-blocking)
    for col in MINIMAL_RECOMMENDED_COLUMNS:
        if col.lower() not in columns_lower:
            logger.warning(
                f"Recommended column '{col}' is missing. Pipeline will use default parameters."
            )

    return errors


@click.command(
    "checksamplesheet",
    short_help="Validate an SDRF file for quantms pipelines.",
)
@click.option("--exp_design", help="SDRF file to be validated", required=True)
@click.option(
    "--template", "-t",
    help="Schema template for full validation (e.g. ms-proteomics, dia-acquisition)",
    default="ms-proteomics",
)
@click.option(
    "--minimal",
    help="Only validate columns required to run the pipeline (skip organism, metadata, etc.)",
    is_flag=True,
)
@click.option(
    "--use_ols_cache_only",
    help="Use OLS cache for ontology validation instead of the live OLS service",
    is_flag=True,
)
def checksamplesheet(
    exp_design: str,
    template: str = "ms-proteomics",
    minimal: bool = False,
    use_ols_cache_only: bool = False,
):
    """Validate an SDRF file for quantms pipelines."""
    check_sdrf(
        input_sdrf=exp_design,
        template=template,
        minimal=minimal,
        use_ols_cache_only=use_ols_cache_only,
    )
