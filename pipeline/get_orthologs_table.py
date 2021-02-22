#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd

DEFAULT_PROJECT_NAME = os.getcwd().split('/')[-1]
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file
Species named as a number
"""


def main(project_name, poff, species, group, genes, required):
    if not poff:
        proteinortho_data = pd.read_csv('{}.proteinortho.tsv'.format(project_name), sep='\t')
        out_file_name = '{}_formed_table.ortho.tsv'.format(project_name)
    else:
        proteinortho_data = pd.read_csv('{}.poff.tsv'.format(project_name), sep='\t')
        out_file_name = '{}_formed_table.ortho.poff.tsv'.format(project_name)


    logging.info("input file {}".format(proteinortho_data))
    species = int(species)
    group = int(group)
    # genes = int(genes)
    # required: to do for multiple required species
    result = proteinortho_data.loc[proteinortho_data['# Species'].isin(range(group, species + 1)) &
                                   (proteinortho_data['{}.faa'.format(required)]!='*')]
    result = result[~result["1.faa"].str.contains(',')]  # because of genes == 1
    result = result[~result["2.faa"].str.contains(',')]  # to do: to normal view
    result = result[~result["3.faa"].str.contains(',')]
    result = result[~result["4.faa"].str.contains(',')]
    result = result[~result["5.faa"].str.contains(',')]
    result = result[~result["6.faa"].str.contains(',')]
    result = result[~result["7.faa"].str.contains(',')]
    result = result[~result["8.faa"].str.contains(',')]
    #  result = proteinortho_data.loc[(proteinortho_data['# Species'] == species)
    #                               & (proteinortho_data['Genes'] == species)]

    with open(out_file_name, 'w') as write_tsv:
        write_tsv.write(result.to_csv(sep='\t', index=False))
        logging.info("out file {} has been recorded".format(out_file_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', help='Project name', nargs='?',
                        default=DEFAULT_PROJECT_NAME)
    parser.add_argument('--poff', help='"y" if use synteny, otherwise empty', nargs='?',
                        default="")
    parser.add_argument('--species', help='Number of species', nargs='?')
    parser.add_argument('--group', help='Minimal size of species group', nargs='?')
    parser.add_argument('--genes', help='Maximal number of genes for species', nargs='?')
    parser.add_argument('--required', help='Required species comma separated', nargs='?')
    args = parser.parse_args()
    try:
        main(args.project, args.poff, args.species, args.group, args.genes, args.required)
    except:
        logging.exception("Unexpected error")

    logging.info("The orthologs table was recorded")
