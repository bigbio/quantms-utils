import logging
import sys

import click
from sdrf_pipelines.sdrf.sdrf import read_sdrf

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_sdrf(
    input_sdrf: str,
    skip_sdrf_validation: bool = False,
    use_ols_cache_only: bool = False,
):
    """
    Check the SDRF file for errors. If any errors are found, print them and exit with a non-zero status code.

    :param input_sdrf: Path to the SDRF file to check
    :param skip_sdrf_validation: Skip all SDRF validation
    :param use_ols_cache_only: Use OLS cache instead of live OLS service
    """
    if skip_sdrf_validation:
        print("No SDRF validation was performed.")
        sys.exit(0)

    df = read_sdrf(input_sdrf)
    errors = df.validate_sdrf(
        template="ms-proteomics",
        use_ols_cache_only=use_ols_cache_only,
    )

    for error in errors:
        print(error)

    sys.exit(bool(errors))



@click.command(
    "checksamplesheet",
    short_help="Validate an SDRF file for quantms pipelines.",
)
@click.option("--exp_design", help="SDRF file to be validated", required=True)
@click.option("--skip_sdrf_validation", help="Skip all SDRF validation", is_flag=True)
@click.option(
    "--use_ols_cache_only",
    help="Use OLS cache for ontology validation instead of the live OLS service",
    is_flag=True,
)
def checksamplesheet(
    exp_design: str,
    skip_sdrf_validation: bool = False,
    use_ols_cache_only: bool = False,
):
    """Validate an SDRF file for quantms pipelines."""
    check_sdrf(
        input_sdrf=exp_design,
        skip_sdrf_validation=skip_sdrf_validation,
        use_ols_cache_only=use_ols_cache_only,
    )
