import argparse
import logging
import os
import re

"""
check duplicates in log files
"""
LOG_FILE = "check_duplicates.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def main(in_dir):
    for infile in os.scandir(in_dir):
        if infile.name.split('.')[-1] == 'log':
            species_names = list()
            with open(os.path.join(in_dir, infile.name), 'r') as f:
                for line in f:
                    if re.search(r'-\s(\d+)$', line):
                        species_names.append(re.search(r'-\s(\d+)$', line).group(1))
            if len(species_names) == len(set(species_names)):
                logging.info("infile {} - OK".format(infile.name))
            else:
                logging.warning("duplicates in file {}".format(infile.name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with .log files to analyze', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder)
    except BaseException as err:
        logging.info("Unexpected error: {}".format(err))
