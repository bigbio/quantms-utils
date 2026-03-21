import logging
import sys

import click

from sdrf_pipelines.sdrf.sdrf import read_sdrf

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def check_sdrf(
    input_sdrf: str,
    template: str = "ms-proteomics",
    skip_sdrf_validation: bool = False,
    use_ols_cache_only: bool = False,
):
    """
    Check the SDRF file for errors. If any errors are found, print them and exit with a non-zero status code.

    :param input_sdrf: Path to the SDRF file to check
    :param template: Schema template to validate against (e.g. 'ms-proteomics', 'dia-acquisition')
    :param skip_sdrf_validation: Skip all SDRF validation
    :param use_ols_cache_only: Use OLS cache instead of live OLS service
    """
    if skip_sdrf_validation:
        print("No SDRF validation was performed.")
        sys.exit(0)

    df = read_sdrf(input_sdrf)
    errors = df.validate_sdrf(
        template=template,
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
@click.option(
    "--template", "-t",
    help="Schema template to validate against (e.g. ms-proteomics, dia-acquisition)",
    default="ms-proteomics",
)
@click.option("--skip_sdrf_validation", help="Skip all SDRF validation", is_flag=True)
@click.option(
    "--use_ols_cache_only",
    help="Use OLS cache for ontology validation instead of the live OLS service",
    is_flag=True,
)
def checksamplesheet(
    exp_design: str,
    template: str = "ms-proteomics",
    skip_sdrf_validation: bool = False,
    use_ols_cache_only: bool = False,
):
    """Validate an SDRF file for quantms pipelines."""
    check_sdrf(
        input_sdrf=exp_design,
        template=template,
        skip_sdrf_validation=skip_sdrf_validation,
        use_ols_cache_only=use_ols_cache_only,
    )
