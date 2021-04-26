#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import traceback
import pandas as pd

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file
Species named as a number
"""


def main(proteinortho_file, species, group, required):
    out_file_name = "{}{}{}".format(os.path.split(proteinortho_file)[0], "/formed_", os.path.split(proteinortho_file)[1])
    proteinortho_data = pd.read_csv(proteinortho_file, sep='\t')
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
    parser.add_argument('--ortho', help='Full path to the output file of Proteinortho', nargs='?')
    parser.add_argument('--species', help='Number of species', nargs='?')
    parser.add_argument('--group', help='Minimal size of species group', nargs='?')
    # parser.add_argument('--genes', help='Maximal number of genes for species', nargs='?') # TODO: is it necessary?
    parser.add_argument('--required', help='One required, target species in relation to which the analysis is made',
                        nargs='?')  # TODO: multi required species?
    args = parser.parse_args()
    try:
        main(args.project, args.ortho, args.species, args.group, args.required)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
    logging.info("The orthologs table was recorded")
