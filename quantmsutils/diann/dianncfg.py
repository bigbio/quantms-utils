"""
This script converts SDRF parameters to DIA-NN parameters
License: Apache 2.0
Authors: Dai Chengxin, Yasset Perez-Riverol
"""

import logging
import re
from typing import List, Tuple
from collections import defaultdict
import click
from sdrf_pipelines.openms.unimod import UnimodDatabase

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)
unimod_database = UnimodDatabase()

@click.command("dianncfg", short_help="Create DIA-NN config file with enzyme and PTMs")
@click.option("--enzyme", "-e", help="")
@click.option("--fix_mod", "-f", help="")
@click.option("--var_mod", "-v", help="")
@click.pass_context
def dianncfg(ctx, enzyme, fix_mod, var_mod):
    """
    Create DIA-NN config file with enzyme and PTMs. It uses the Unimod database to validate and format the
    modifications for DIA-NN compatibility. These strings, along with the enzyme cut rule, are written to a
    configuration file named 'diann_config.cfg'.

    :param ctx: Click context
    :param enzyme: The name of the enzyme used for protein digestion.
    :param fix_mod: A string of fixed modifications, separated by commas.
    :param var_mod: A string of variable modifications, separated by commas.
    """
    cut = enzyme_cut(enzyme)
    fix_ptm, var_ptm = convert_mod(fix_mod, var_mod)

    var_ptm_str = " --var-mod "
    fix_ptm_str = " --fixed-mod "
    diann_fix_ptm = ""
    diann_var_ptm = ""
    for mod in fix_ptm:
        diann_fix_ptm += fix_ptm_str + mod
    for mod in var_ptm:
        diann_var_ptm += var_ptm_str + mod

    with open("diann_config.cfg", "w") as file:
        file.write("--cut " + cut + diann_fix_ptm + diann_var_ptm)


def get_mod(mod):
    pattern = re.compile(r"\((.*?)\)")
    tag = 0
    diann_mod = None
    for modification in unimod_database.modifications:
        if modification.get_name() == mod.split(" ")[0]:
            diann_mod = modification.get_name() + "," + str(modification._delta_mono_mass)
            tag = 1
            break
    if tag == 0:
        logging.error(
            "Currently only supported unimod modifications for DIA pipeline. Skipped: "
            + mod
        )
    site = re.findall(pattern, " ".join(mod.split(" ")[1:]))[0]
    if site == "Protein N-term":
        site = "*n"
    elif site == "N-term":
        site = "n"

    if (
            "TMT" in diann_mod
            or "Label:" in diann_mod
            or "iTRAQ" in diann_mod
            or "mTRAQ" in diann_mod
            or "Dimethyl:" in diann_mod
    ):
        logging.error(
            "quantms DIA-NN workflow only support LFQ now!"
            + mod
        )
    elif diann_mod is not None:
        return diann_mod, site
    else:
        logging.error(
            "Currently only supported unimod modifications for DIA pipeline. Skipped: "
            + mod
        )
        return


def convert_mod(fix_mod: str, var_mod: str) -> Tuple[List, List]:
    var_ptm = []
    fix_ptm = []
    if fix_mod != "":
        merged = defaultdict(list)
        for mod in fix_mod.split(","):
            diann_mod, site = get_mod(mod)
            merged[diann_mod].append(site)

        # merge same modification for different sites
        for name, site_list in merged.items():
            site_str = "".join(sorted(set(site_list)))
            fix_ptm.append(f"{name} ({site_str})")

    if var_mod != "":
        merged = defaultdict(list)
        for mod in var_mod.split(","):
            diann_mod, site = get_mod(mod)
            merged[diann_mod].append(site)
        # merge same modification for different sites
        for name, site_list in merged.items():
            site_str = "".join(sorted(set(site_list)))
            var_ptm.append(f"{name} ({site_str})")

    return fix_ptm, var_ptm


_ENZYME_SPECIFICITY = {
    "Trypsin": "K*,R*,!*P",
    "Trypsin/P": "K*,R*",
    "Arg-C": "R*,!*P",
    "Asp-N": "*B,*D",
    "Chymotrypsin": "F*,W*,Y*,L*,!*P",
    "Lys-C": "K*,!*P",
}


def enzyme_cut(enzyme: str) -> str:
    return _ENZYME_SPECIFICITY.get(enzyme) or "--cut"
