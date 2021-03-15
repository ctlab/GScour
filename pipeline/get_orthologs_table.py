#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import traceback

import pandas as pd

DEFAULT_PROJECT_NAME = os.getcwd().split('/')[-1]
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file
Species named as a number
"""


def main(project_name, poff, species, group, required):
    if not poff:
        proteinortho_file = '{}.proteinortho.tsv'.format(project_name)
        proteinortho_data = pd.read_csv(proteinortho_file, sep='\t')
        out_file_name = '{}_formed_table.ortho.tsv'.format(project_name)
    else:
        proteinortho_file = '{}.poff.tsv'.format(project_name)
        proteinortho_data = pd.read_csv(proteinortho_file, sep='\t')
        out_file_name = '{}_formed_table.ortho.poff.tsv'.format(project_name)

    logging.info("Input file {}".format(proteinortho_file))
    species = int(species)
    group = int(group)
    result = proteinortho_data.loc[proteinortho_data['# Species'].isin(range(group, species + 1)) &
                                   (proteinortho_data['{}.faa'.format(required)] != '*')]
    for species_column in result.columns[3: 3 + species + 1]:
        result = result[~result[species_column].str.contains(',')]  # because of genes == 1 (single-copy orthologs)

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
    # parser.add_argument('--genes', help='Maximal number of genes for species', nargs='?') # TODO: is it necessary?
    parser.add_argument('--required', help='One required, target species', nargs='?')  # TODO: multi required species?
    args = parser.parse_args()
    try:
        main(args.project, args.poff, args.species, args.group, args.required)
    except BaseException as e:
        logging.info("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))

    logging.info("The orthologs table was recorded")
