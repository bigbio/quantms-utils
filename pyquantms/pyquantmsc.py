import click

from pyquantms.diann.create_diann_cfg import create_diann_cfg

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(create_diann_cfg)

if __name__ == "__main__":
    cli()
