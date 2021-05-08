#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
import numpy as np

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file in accordance with predefined list of genes.
Predefined list of genes wrote in .xlsx with one column = 'Gene name'.
Species named as a number.
"""


def get_protein_ids_from_gbff(gb_file, gene_names, protein_length_dict, protein_genes_dict):
    """parse annotation file for target species to get protein_ids by gene_name"""
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("gene") in gene_names.values:
                    protein_id = feature.qualifiers.get("protein_id")[0]  # str
                    protein_translation = feature.qualifiers.get("translation")[0]  # str
                    protein_translation_length = len(protein_translation)  # int
                    gene_name = feature.qualifiers.get('gene')[0]
                    protein_genes_dict[protein_id] = gene_name
                    protein_info = [{'protein_id': protein_id, 'length': protein_translation_length}]
                    if gene_name not in protein_length_dict.keys():
                        protein_length_dict[gene_name] = list()
                    protein_length_dict[gene_name] += protein_info


def write_to_ortho_file(ortho_file_path, protein_genes_dict, required_species):
    ortho_data = pd.read_csv(ortho_file_path, sep='\t')
    required_column = '{}.{}'.format(required_species, 'faa')
    gene_names_column = []
    for protein_id in ortho_data[required_column]:
        if protein_id in protein_genes_dict.keys():
            gene_names_column.append(protein_genes_dict[protein_id])
        else:
            gene_names_column.append(np.NaN)
    ortho_data['Gene_name'] = gene_names_column
    ortho_data.to_csv(ortho_file_path, sep='\t')


def main(ortho_file, gb_file_path, predefined_genes_path, required_species):
    log_file = os.path.join(os.path.split(ortho_file)[0], "{}.{}".format('predefined', 'log'))
    child_logger = logging.getLogger(__name__)
    child_logger.addHandler(logging.FileHandler(log_file))
    child_logger.setLevel(10)
    child_logger.info("Work with args:\n ortho_file {}\n gen bank annotation file {}\n"
                      "predefined genes file {}\n required species number {}\n"
                      "log file will be written to {}".format(ortho_file, gb_file_path, predefined_genes_path,
                                                              required_species, log_file))

    protein_length_dict = dict()
    protein_genes_dict = dict()
    df = pd.io.excel.read_excel(predefined_genes_path)
    get_protein_ids_from_gbff(gb_file_path, df['Gene name'], protein_length_dict, protein_genes_dict)
    child_logger.info("Common dict with proteins length of length {}:\n{}".
                      format(len(protein_length_dict), protein_length_dict))
    child_logger.info("Protein-gene dict to csv write of length {}:\n{}".
                      format(len(protein_genes_dict), protein_genes_dict))
    write_to_ortho_file(ortho_file, protein_genes_dict, required_species)
    child_logger.info("Genes added to ortho data file. Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Full path to the output file of some ortho finder (will be modified,'
                                        'do a copy is recommended)', nargs='?')
    parser.add_argument('--gbff', help='Path to the .gbff annotation file for target species', nargs='?')
    parser.add_argument('--predefined', help='Path to the .xlsx file with one 1st column - \'Gene symbol\'',
                        nargs='?')
    parser.add_argument('--required', help='One required, number of target species in relation to which the analysis'
                                           ' is made', nargs='?')
    args = parser.parse_args()
    try:
        main(args.ortho, args.gbff, args.predefined, args.required)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))

