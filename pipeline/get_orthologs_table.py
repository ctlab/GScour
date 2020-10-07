import argparse
import logging
import os
import pandas as pd

DEFAULT_PROJECT_NAME = os.getcwd().split('/')[-1]


def main(project_name, species):
    proteinortho_data = pd.read_csv('{}.proteinortho.tsv'.format(project_name), sep='\t')
    species = int(species)
    result = proteinortho_data.loc[(proteinortho_data['# Species'] == species)
                                   & (proteinortho_data['Genes'] == species)]
    with open('single_copy_orthologs.tsv','w') as write_tsv:
        write_tsv.write(result.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', help='Project name', nargs='?',
                        default=DEFAULT_PROJECT_NAME)
    parser.add_argument('--species', help='Number of species', nargs='?')
    args = parser.parse_args()
    try:
        main(args.project, args.species)
    except:
        logging.exception("Unexpected error")

    logging.info("The orthologs table was recorded")
