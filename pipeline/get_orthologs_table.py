#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file
Species named as a number
"""


def get_gene_names(df, gbff_path, required_species):
    """parse annotation file to get gene_name by protein_id"""
    protein_gene = list()
    for record in SeqIO.parse(gbff_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("protein_id") in df[required_species].values:
                    gene_name = feature.qualifiers.get('gene')[0]
                    protein_id = feature.qualifiers.get("protein_id")[0]
                    protein_gene.append([protein_id, gene_name])
    return pd.DataFrame(protein_gene, columns=[required_species, 'Gene name'])


def main(proteinortho_file, species, required, gbff_path):
    out_file_name = "{}{}{}".format(os.path.split(proteinortho_file)[0], "/formed_", os.path.split(proteinortho_file)[1])
    proteinortho_data = pd.read_csv(proteinortho_file, sep='\t')
    logging.info("Input file {}\nOutput file {}".format(proteinortho_file, out_file_name))
    result = []
    for species in species.split(','):
        species = int(species)
        res = proteinortho_data.loc[proteinortho_data['# Species'].isin(range(species, species + 1)) &
                                    proteinortho_data['Genes'].isin(range(species, species + 1))]
        result.append(res)
    result = pd.concat(result)
    result = result.loc[result['{}'.format(required)] != '*']
    prot_gene_df = get_gene_names(result, gbff_path, required)
    result = result.merge(prot_gene_df, on=required)
    with open(out_file_name, 'w') as f:
        f.write(result.to_csv(sep='\t'))
        logging.info("Output file {} has been recorded".format(out_file_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Full path to the output file of Proteinortho', nargs='?', required=True)
    parser.add_argument('--species', help='Int values for column \'# Species\'(\'Genes\') separated by comma'
                                          'relationship between species and genes, because of searching one-to-one '
                                          'orthologs. For example \'3, 4, 5\'', nargs='?', required=True)
    parser.add_argument('--required', help='One required, number of target species in relation to which the analysis'
                                           ' is made', nargs='?', required=True)  # TODO: multi required species?
    parser.add_argument('--gbff', help='Path to the .gbff annotation file to extract gene names', nargs='?',
                        required=True)
    args = parser.parse_args()
    try:
        main(args.ortho, args.species, args.required, args.gbff)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    logging.info("The orthologs table was recorded")
