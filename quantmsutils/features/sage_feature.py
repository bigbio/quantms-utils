import click
import pandas as pd
import pyopenms as oms


@click.command("sage2feature")
@click.option("--idx_file", "-i", help="Input idXML file")
@click.option("--output_file", "-o", help="Output idXML file")
@click.option("--feat_file", "-f", help="Input feature table file")
@click.pass_context
def add_sage_feature(ctx, idx_file: str, output_file: str, feat_file: str):
    """
    Add extra features in features idXML. Adding extra feature in Sage isn't known input for PSMFeatureExtractor

    :param ctx: click context
    :param idx_file: Original idXML file
    :param output_file: Outpuf file with the extra feature
    :param feat_file: Feature file from Sage
    :return: None
    """
    extra_feat = []
    feat = pd.read_csv(feat_file, sep="\t")
    for _, row in feat.iterrows():
        if row["feature_generator"] == "psm_file":
            continue
        else:
            extra_feat.append(row["feature_name"])
    print("Adding extra feature: {}".format(extra_feat))
    protein_ids = []
    peptide_ids = []
    oms.IdXMLFile().load(idx_file, protein_ids, peptide_ids)
    search_parameters = protein_ids[0].getSearchParameters()
    features = search_parameters.getMetaValue("extra_features")
    extra_features = features + "," + ",".join(extra_feat)
    search_parameters.setMetaValue("extra_features", extra_features)
    protein_ids[0].setSearchParameters(search_parameters)
    oms.IdXMLFile().store(output_file, protein_ids, peptide_ids)
    print("Done")
