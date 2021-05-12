#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
import operator

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file 
Reduce redundant number of proteins to single one (choose the longest)
Species named as a number
"""


def get_longest_protein(gbff_path, column, list_of_proteins, gene_name):
    """ parse annotation file to get protein_id of the longest protein
        column name e.g. '6.faa', 6 - number of species
        6.gbff name of genbank annotation file
    """
    protein_length_dict = dict()
    for record in SeqIO.parse(os.path.join(gbff_path, '{}.{}'.format(column.split('.')[0], 'gbff')), "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and feature.qualifiers.get("protein_id"):
                protein_id = feature.qualifiers.get("protein_id")[0]  # str
                if protein_id in list_of_proteins:
                    protein_translation = feature.qualifiers.get("translation")[0]  # str
                    protein_translation_length = len(protein_translation)  # int
                    protein_length_dict[protein_id] = protein_translation_length
                    gene_name = feature.qualifiers.get("gene")
    if protein_length_dict:
        logging.info("protein_length_dict for gene {}:{}".format(gene_name, protein_length_dict))
        longest_protein_id = max(protein_length_dict.items(), key=operator.itemgetter(1))[0]
        return longest_protein_id
    else:
        logging.info("nothing find for gene {}, proteins {}".format(gene_name, list_of_proteins))
        return


def main(proteinortho_file, gbff_path):
    proteinortho_data = pd.read_csv(proteinortho_file, sep='\t')
    logging.info("Input file {}".format(proteinortho_file))
    for idx, row in proteinortho_data.iterrows():
        for proteins in row:
            for column in row.keys():
                if row[column] == proteins:
                    target_column = column
                    gene_name = row['Gene_name']
            if isinstance(proteins, str) and ',' in proteins:
                list_of_proteins = proteins.split(',')
                longest_protein = get_longest_protein(gbff_path, target_column, list_of_proteins, gene_name)
                if longest_protein:
                    proteinortho_data[target_column] = proteinortho_data[target_column].replace([proteins],
                                                                                                longest_protein)
                    proteinortho_data['Genes'].iloc[idx] = proteinortho_data['Genes'].iloc[idx] \
                                                           - len(list_of_proteins) - 1
            # elif isinstance(proteins, str):  # just for information
            #     protein = [proteins]
            #     get_longest_protein(gbff_path, target_column, protein, gene_name)

    with open(proteinortho_file, 'w') as write_tsv:
        write_tsv.write(proteinortho_data.to_csv(sep='\t', index=False))
        logging.info("The input file is overwritten")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Full path to the output file of Proteinortho, will be overwritten do the copy'
                                        'is recommended', nargs='?')
    parser.add_argument('--gbff', help='Path to the folder with .gbff annotation files for every species', nargs='?')
    args = parser.parse_args()
    try:
        main(args.ortho, args.gbff)  # TODO: multiprocessing
    except BaseException as e:
        logging.info("Unexpected error: {}".format(e))
        raise e
    logging.info("The orthologs table with many-to-many orthologs reduced to one-to-one orthologs.")
