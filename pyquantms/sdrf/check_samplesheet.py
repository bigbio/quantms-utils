# nf-core: Update the script to check the sdrf
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import argparse
import errno
import os
import sys

import pandas as pd
from sdrf_pipelines.sdrf.sdrf import SdrfDataFrame
from sdrf_pipelines.sdrf.sdrf_schema import DEFAULT_TEMPLATE, MASS_SPECTROMETRY


def parse_args(args=None):
    description = "Reformat nf-core/quantms sdrf file and check its contents."
    epilog = "Example usage: python validate_sdrf.py <sdrf> <check_ms>"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("SDRF", help="SDRF/Expdesign file to be validated")
    parser.add_argument("ISSDRF", help="SDRF file or Expdesign file")
    parser.add_argument(
        "--CHECK_MS",
        help="check mass spectrometry fields in SDRF.",
        action="store_true",
    )

    return parser.parse_args(args)


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


def check_sdrf(check_ms, sdrf):
    df = SdrfDataFrame.parse(sdrf)
    errors = df.validate(DEFAULT_TEMPLATE)
    if check_ms:
        errors = errors + df.validate(MASS_SPECTROMETRY)
    for error in errors:
        print(error)
    if not errors:
        print("Everying seems to be fine. Well done.")
    else:
        print("There were validation errors!")
    sys.exit(bool(errors))


def check_expdesign(expdesign):
    data = pd.read_csv(expdesign, sep="\t", header=0, dtype=str)
    data = data.dropna()
    schema_file = ["Fraction_Group", "Fraction", "Spectra_Filepath", "Label", "Sample"]
    schema_sample = ["Sample", "MSstats_Condition", "MSstats_BioReplicate"]

    # check table format: two table
    with open(expdesign, "r") as f:
        lines = f.readlines()
        try:
            empty_row = lines.index("\n")
        except ValueError:
            print(
                "the one-table format parser is broken in OpenMS2.5, please use one-table or sdrf"
            )
            sys.exit(1)
        if lines.index("\n") >= len(lines):
            print(
                "the one-table format parser is broken in OpenMS2.5, please use one-table or sdrf"
            )
            sys.exit(1)

        s_table = [i.replace("\n", "").split("\t") for i in lines[empty_row + 1 :]][1:]
        s_header = lines[empty_row + 1].replace("\n", "").split("\t")
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    # check missed mandatory column
    missed_columns = set(schema_file) - set(data.columns)
    if len(missed_columns) != 0:
        print("{0} column missed".format(" ".join(missed_columns)))
        sys.exit(1)

    missed_columns = set(schema_sample) - set(s_data_frame.columns)
    if len(missed_columns) != 0:
        print("{0} column missed".format(" ".join(missed_columns)))
        sys.exit(1)

    if len(set(data.Label)) != 1 and "MSstats_Mixture" not in s_data_frame.columns:
        print("MSstats_Mixture column missed in ISO experiments")
        sys.exit(1)

    # check logical problem: may be improved
    check_expdesign_logic(data, s_data_frame)


def check_expdesign_logic(f_table, s_table):
    if int(max(f_table.Fraction_Group)) > len(set(f_table.Fraction_Group)):
        print("Fraction_Group discontinuous!")
        sys.exit(1)
    f_table_d = f_table.drop_duplicates(
        ["Fraction_Group", "Fraction", "Label", "Sample"]
    )
    if f_table_d.shape[0] < f_table.shape[0]:
        print(
            "Existing duplicate entries in Fraction_Group, Fraction, Label and Sample"
        )
        sys.exit(1)
    if len(set(s_table.Sample)) < s_table.shape[0]:
        print("Existing duplicate Sample in sample table!")
        sys.exit(1)


def main(args=None):
    # TODO validate expdesign file
    args = parse_args(args)

    if args.ISSDRF == "true":
        check_sdrf(args.CHECK_MS, args.SDRF)
    else:
        check_expdesign(args.SDRF)


if __name__ == "__main__":
    sys.exit(main())
