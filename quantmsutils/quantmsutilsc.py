import click

from quantmsutils.diann.dianncfg import dianncfg
from quantmsutils.diann.diann2mztab import diann2mztab
from quantmsutils.features.sage_feature import add_sage_feature
from quantmsutils.rescoring.ms2rescore import ms2rescore
from quantmsutils.mzml.mzml_statistics import mzml_statistics
from quantmsutils.psm.psm_conversion import convert_psm
from quantmsutils.sdrf.check_samplesheet import check_samplesheet
from quantmsutils.sdrf.extract_sample import extract_sample_from_expdesign

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(dianncfg)
cli.add_command(diann2mztab)

cli.add_command(add_sage_feature)

cli.add_command(mzml_statistics)

cli.add_command(extract_sample_from_expdesign)
cli.add_command(check_samplesheet)

cli.add_command(ms2rescore)

cli.add_command(convert_psm)

def main():
    try:
        cli()
    except SystemExit as e:
        if e.code != 0:
            raise
if __name__ == "__main__":
    main()
