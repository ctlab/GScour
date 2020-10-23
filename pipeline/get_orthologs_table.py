import argparse
import logging
import os
import pandas as pd

DEFAULT_PROJECT_NAME = os.getcwd().split('/')[-1]


def main(project_name, species, poff):
    if not poff:
        proteinortho_data = pd.read_csv('{}.proteinortho.tsv'.format(project_name), sep='\t')
        out_file_name = 'single_copy_orthologs_{}.tsv'.format(project_name)
    else:
        proteinortho_data = pd.read_csv('{}.poff.tsv'.format(project_name), sep='\t')
        out_file_name = 'single_copy_orthologs_poff_{}.tsv'.format(project_name)
    species = int(species)
    result = proteinortho_data.loc[(proteinortho_data['# Species'] == species)
                                   & (proteinortho_data['Genes'] == species)]
    with open(out_file_name, 'w') as write_tsv:
        write_tsv.write(result.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', help='Project name', nargs='?',
                        default=DEFAULT_PROJECT_NAME)
    parser.add_argument('--poff', help='"y" if use synteny, otherwise empty', nargs='?',
                        default="")
    parser.add_argument('--species', help='Number of species', nargs='?')
    args = parser.parse_args()
    try:
        main(args.project, args.species, args.poff)
    except:
        logging.exception("Unexpected error")

    logging.info("The orthologs table was recorded")
