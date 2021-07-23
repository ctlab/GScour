#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
"""
Form ortholog table from proteinortho.tsv file
Species named as a number
"""


def main(proteinortho_file, species, required):
    out_file_name = "{}{}{}".format(os.path.split(proteinortho_file)[0], "/formed_", os.path.split(proteinortho_file)[1])
    proteinortho_data = pd.read_csv(proteinortho_file, sep=',')
    logging.info("Input file {}\nOutput file {}".format(proteinortho_file, out_file_name))
    result = []
    for species in species.split(','):
        species = int(species)
        res = proteinortho_data.loc[proteinortho_data['# Species'].isin(range(species, species + 1)) &
                                    proteinortho_data['Genes'].isin(range(species, species + 1))]
        result.append(res)
    result = pd.concat(result)
    result = result.loc[result['{}.faa'.format(required)] != '*']
    # result = proteinortho_data.loc[proteinortho_data['# Species'].isin(range(group, species + 1)) &
    #                                (proteinortho_data['{}.faa'.format(required)] != '*')]
    # for species_column in result.columns[3: 3 + species]:
    #     result = result[~result[species_column].str.contains(',')]  # because of genes == 1 (single-copy orthologs)

    #  result = proteinortho_data.loc[(proteinortho_data['# Species'] == species)
    #                               & (proteinortho_data['Genes'] == species)]

    with open(out_file_name, 'w') as write_tsv:
        write_tsv.write(result.to_csv(sep='\t', index=False))
        logging.info("Output file {} has been recorded".format(out_file_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Full path to the output file of Proteinortho', nargs='?', required=True)
    parser.add_argument('--species', help='Int values for column \'# Species\' separated by comma are equal '
                                          'Int values for column \'Genes\', assuming one-to-one '
                                          'relationship between species and genes, because of searching one-to-one '
                                          'orthologs', nargs='?', required=True)
    # parser.add_argument('--genes', help='Int values for column \'Genes\' separated by comma', nargs='?')
    parser.add_argument('--required', help='One required, number of target species in relation to which the analysis'
                                           ' is made', nargs='?', required=True)  # TODO: multi required species?
    args = parser.parse_args()
    try:
        main(args.ortho, args.species, args.required)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    logging.info("The orthologs table was recorded")
