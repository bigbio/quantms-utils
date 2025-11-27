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

# Lazy initialization of UnimodDatabase for improved testability.
# The database is created on first access rather than at module import time,
# which allows tests to mock or replace it more easily.
_unimod_database = None


def get_unimod_database():
    """
    Get the UnimodDatabase instance, creating it lazily on first access.

    This pattern improves testability by avoiding database initialization at module
    import time. For testing purposes, the internal _unimod_database variable can be
    set to None to force re-initialization on the next call.

    :return: The UnimodDatabase instance.
    """
    global _unimod_database
    if _unimod_database is None:
        _unimod_database = UnimodDatabase()
    return _unimod_database


# Met-loss modification constant (UniMod:765) with mass shift and site specification
MET_LOSS_MODIFICATION = "UniMod:765,-131.040485,*nM"


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
        if mod == MET_LOSS_MODIFICATION:
            diann_var_ptm += " --met-excision "
        else:
            diann_var_ptm += var_ptm_str + mod

    with open("diann_config.cfg", "w") as file:
        file.write("--cut " + cut + diann_fix_ptm + diann_var_ptm)


def get_mod(mod, mod_type):
    """
    Retrieve and format a modification from the Unimod database for DIA-NN compatibility.

    :param mod: The modification string, typically containing the modification name and site.
    :param mod_type: The type of modification ('fixed_mod' or 'var_mod').
    :return: A tuple (diann_mod_accession, site), where diann_mod_accession is a formatted string
             for DIA-NN and site is the modification site.
    :raises SystemExit: If the modification is not found in the Unimod database, logs an error and exits.
    """
    pattern = re.compile(r"\((.*?)\)")
    modification_found = 0
    diann_mod_accession = None
    diann_mod_name = None
    for modification in get_unimod_database().modifications:
        if modification.get_name() == mod.split(" ")[0]:
            diann_mod_accession = modification.get_accession().replace("UNIMOD:", "UniMod:") + "," + str(modification._delta_mono_mass)
            diann_mod_name = modification.get_name()
            modification_found = 1
            break

    if modification_found == 0:
        logging.error(
            f"Only Unimod modifications are currently supported for the DIA pipeline. Unsupported modification: {mod}"
        )
        exit(1)

    # TODO support DIA multiplex
    if (
            "TMT" in diann_mod_name
            or "Label:" in diann_mod_name
            or "iTRAQ" in diann_mod_name
            or "mTRAQ" in diann_mod_name
            or "Dimethyl:" in diann_mod_name
    ):
        logging.error(
            "quantms DIA-NN workflow only supports LFQ now! Unsupported modifications: "
            + mod
        )
        exit(1)

    sites = re.findall(pattern, " ".join(mod.split(" ")[1:]))
    if not sites:
        logging.error(
            f"No site specification found in modification string: {mod}"
        )
        exit(1)
    site = sites[0]
    if site == "Protein N-term":
        site = "*n"
    elif site == "N-term":
        site = "n"
    elif len(site.split(" ")) >= 2:
        pp = " ".join(site.split(" ")[:-1])
        if pp == "Protein N-term":
            pp = "*n"
        elif pp == "N-term":
            pp = "n"
        aa = site.split(" ")[-1]
        site = pp + aa
        if site == "*nM" and diann_mod_name == "Met-loss" and mod_type == "var_mod":
            return diann_mod_accession, site
        else:
            logging.warning("Restricting to certain terminal AAs isn't directly supported. Please see https://github.com/vdemichev/DiaNN/issues/1791")
    return diann_mod_accession, site


def convert_mod(fix_mod: str, var_mod: str) -> Tuple[List, List]:
    var_ptm = []
    fix_ptm = []
    if fix_mod:
        merged = defaultdict(list)
        for mod in fix_mod.split(","):
            diann_mod, site = get_mod(mod, "fixed_mod")
            merged[diann_mod].append(site)

        # merge same modification for different sites
        for name, site_list in merged.items():
            site_str = "".join(sorted(set(site_list)))
            fix_ptm.append(f"{name},{site_str}")

    if var_mod:
        merged = defaultdict(list)
        for mod in var_mod.split(","):
            diann_mod, site = get_mod(mod, "var_mod")
            merged[diann_mod].append(site)
        # merge same modification for different sites
        for name, site_list in merged.items():
            site_str = "".join(sorted(set(site_list)))
            var_ptm.append(f"{name},{site_str}")

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
