#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os

"""remove broken items
broken items wrote down to dict BROKEN_LIST
BROKEN_LIST = {'item_folder':[item1, item2]}
"""
BROKEN_LIST = {}
NOT_FOUND_FILES = list()


def replace_broken_files(directory_out):
    global NOT_FOUND_FILES
    broken_multiple_folder = "broken_multiple_files"
    if not os.path.isdir(directory_out):
        os.makedirs(broken_multiple_folder)
    try:
        for file_number in BROKEN_LIST.keys():
            os.replace(os.path.join(directory_out, file_number + ".fna"), os.path.join(broken_multiple_folder, file_number + ".fna"))
            os.replace(os.path.join(directory_out, file_number + ".log"), os.path.join(broken_multiple_folder, file_number + ".log"))
            logging.info("Replaced all from BROKEN_LIST to the folder {}".format(broken_multiple_folder))
    except FileNotFoundError:
        if file_number not in NOT_FOUND_FILES:
            NOT_FOUND_FILES.append(file_number)


def main(outfolder):
    replace_broken_files(outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfolder', help='Path to the folder with broken files', nargs='?')
    args = parser.parse_args()
    try:
        main(args.outfolder)
    except BaseException as err:
        logging.info("Unexpected error: {}".format(err))
    logging.info("NOT_FOUND_FILES list of length {}:\n{}".format(len(NOT_FOUND_FILES), NOT_FOUND_FILES))


