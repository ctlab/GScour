#!/usr/bin/env python
import argparse
import os
import logging
import traceback

"""
script for launch swamp (python3)
target_dict - if there is certain files to launch
target_dict in format:
target_dict[species_folder] = [item_folder1, item_folder2...]
"""
LOG_FILE = os.path.join(os.getcwd(), 'swamp_log.log')
print("LOG_FILE", LOG_FILE)
BROKEN_FILES = list()
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
# TODO: child_logger
target_dict = {}


def get_branch_names_file_path(branch_names_folder, species_folder_name):
    for branch_name_file in os.scandir(branch_names_folder):
        if branch_name_file.name.split('.')[0] == species_folder_name:
            return branch_name_file.name


def get_input_items(folder_in, branch_name_folder):
    """ parse root folder with files for paml
    parse branch_names_folder to get appropriate branches code file
    return item folder and path to the branch_code file """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):
            branch_name = get_branch_names_file_path(branch_name_folder, species_folder.name)
            branch_name_path = os.path.join(branch_name_folder, branch_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if infile.split('.')[-1] == 'phy':
                            yield os.path.join(folder_in, species_folder.name, item.name), branch_name_path


def get_target_input_items(folder_in, branch_name_folder):
    """ parse root folder with files for paml
    parse branch_names_folder to get appropriate branch_code file
    return item folder and path to the branch_code file
    return item folder if it's in the target_dict """
    global target_dict
    for species_folder in os.scandir(folder_in):
        species_folder_path = os.path.join(folder_in, species_folder.name)
        if os.path.isdir(species_folder_path):
            if target_dict.get(species_folder.name):
                branch_name = get_branch_names_file_path(branch_name_folder, species_folder.name)
                branch_name_path = os.path.join(branch_name_folder, branch_name)
                for item in os.scandir(species_folder_path):
                    item_folder_path = os.path.join(folder_in, species_folder.name, item.name)
                    if os.path.isdir(item_folder_path):
                        if item.name in target_dict[species_folder.name]:
                            for in_file in os.scandir(item_folder_path):
                                if in_file.name.split('.')[-1] == 'phy':
                                    yield os.path.join(folder_in, species_folder.name, item.name), branch_name_path


def run_swamp(items_folder, executable_path, branch_codes, threshold, windows_size):
    global BROKEN_FILES
    launch_swamp = '{}' \
                   ' -i {} -b {} -t {} -w {} >> {}'.format(executable_path, items_folder, branch_codes,
                                                           threshold, windows_size, LOG_FILE)
    try:
        if os.system(launch_swamp):  # TODO: catching full swamp stderr
            raise ValueError
    except ValueError as err:
        file_number = items_folder.split('/')[-1]
        logging.exception("File {} error args: {}".format(items_folder, err.args))
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
        if target_dict:
            for item_folder, branch_names_file in get_target_input_items(args.i, args.b):
                run_swamp(item_folder, args.e, branch_names_file, args.t, args.w)
        else:
            for item_folder, branch_names_file in get_input_items(args.i, args.b):
                run_swamp(item_folder, args.e, branch_names_file, args.t, args.w)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
    logging.info("BROKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    logging.info("The work has been completed")
