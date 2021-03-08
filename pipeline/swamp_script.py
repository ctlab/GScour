#!/usr/bin/env python
import argparse
import re
import sys
import logging
import os
import logging

LOG_FILE = '../usages/swamp_log.log'
BROKEN_FILES = list()
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def get_branch_names_file_path(branch_names_folder, species_folder_name):
    for branch_name_file in os.scandir(branch_names_folder):
        if branch_name_file.name.split('.')[0] == species_folder_name:
            return branch_name_file.name


def get_input_items(folder_in, branch_name_folder):
    """ parse root folder with files for paml
    parse branch_names_folder to get appropriate branches code file
    return item folder and path to the file with branch codes"""
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):
            branch_name = get_branch_names_file_path(branch_name_folder, species_folder.name)
            branch_name_path = os.path.join(branch_name_folder, branch_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if infile.split('.')[-1] == 'phy':
                            yield os.path.join(folder_in, species_folder.name, item.name), branch_name_path


def run_swamp(items_folder, executable_path, branch_codes, threshold, windows_size):
    global BROKEN_FILES
    try:
        launch_swamp = '{}' \
                       ' -i {} -b {} -t {} -w {} >> {}'.format(executable_path, items_folder, branch_codes,
                                                               threshold, windows_size, LOG_FILE)
        if os.system(launch_swamp):
            raise ValueError      # TODO: catch the full stderr from swamp; target dict
    except ValueError:
        logging.exception("sys.exc_info() {0}".format(sys.exc_info()))
        file_number = (re.search(r"/(\d+)/*", items_folder)).group(1)
        if file_number not in BROKEN_FILES:
            BROKEN_FILES.append(file_number)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', help='Path to the SWAMP executable', nargs='?')
    parser.add_argument('-i', help='Provide the full path to an INFOLDER that contain multiple subfolders',
                        nargs='?')
    parser.add_argument('-b', help='Path to the folder with BRANCHNAMESFILEs', nargs='?')
    parser.add_argument('-t', help='THRESHOLD', nargs='?', default="")
    parser.add_argument('-w', help='WINDOWSIZE', nargs='?')
    args = parser.parse_args()
    try:
        for item_folder, branch_names_file in get_input_items(args.i, args.b):
            run_swamp(item_folder, args.e, branch_names_file, args.t, args.w)
    except:
        logging.exception("Unexpected error")
    logging.info("BROKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    logging.info("The work has been completed")