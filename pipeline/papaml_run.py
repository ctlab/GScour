#!/usr/bin/env python
import argparse
import os
import logging

CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()

LOG_FILE = "papaml.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)

"""
see https://github.com/RetroWWU/paPAML
"""


def get_input_items(folder_in):
    """ parse root folder with files for paml
    parse tree_folder to get appropriate tree """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):  # and species_folder.name.isdigit():
            logging.info("working with species folder {}".format(species_folder.name))
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    yield species_folder.name, item.name


def main(folder_in, exec_path, number_of_threads):
    for species_folder, item in get_input_items(folder_in):
        os.chdir(os.path.join(in_folder, species_folder, item))
        logging.info("working with item {}".format(item))
        os.system('perl {}paPAML.pl -p {} -d'.format(exec_path, number_of_threads))
        # os.system('perl {}paPAML.pl -c'.format(exec_path)) # for cleaning folder


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', help='Path to the papaml executable', nargs='?', default="paPAML.pl")
    # parser.add_argument('--timeout', help='Timeout for papaml in seconds, default=500', nargs='?', default='500')
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?', required=True)
    parser.add_argument('--threads', help='Number of threads to use', nargs='?', required=True)
    # parser.add_argument('--rework', help='"y" if overwrite existing files, default "n"', nargs='?', default='n')
    args = parser.parse_args()
    in_folder = args.i
    executable_path = args.e
    threads = int(args.threads)
    # timeout = int(args.timeout)
    # if args.rework == 'y':
    #     rework = True
    # else:
    #     rework = False
    logging.info("Path to the folder with input files for papaml: {}\nExecutable path: {}\n"
                 "Threads to use = {}".
                 format(in_folder, executable_path, threads))
    try:
        main(in_folder, executable_path, threads)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    logging.info("The work has been completed")
