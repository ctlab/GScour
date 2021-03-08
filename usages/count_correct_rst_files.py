#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os

LOG_FILE = "count_rst.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def get_rst_file(folder_in):
    """ parse root folder with files for paml
    to get rst file """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if infile == 'rst':
                            yield folder_in, species_folder.name, item.name, infile


def main(folder_in):
    right_rst_counter = 0
    rst_dict = dict()
    counter = 0
    counter_empty = 0
    for folder_in, species_folder, item_folder, infile in get_rst_file(folder_in):
        counter += 1
        rst_file = os.path.join(folder_in, species_folder, item_folder, infile)
        with open(rst_file, 'r') as f:
            if os.path.getsize(rst_file) > 0:
                for line in f:
                    if "Summary of changes along branches" in line:
                        right_rst_counter += 1
                        if not rst_dict.get(species_folder):
                            rst_dict[species_folder] = list()
                        rst_dict[species_folder].append(item_folder)
            else:
                counter_empty += 1

    logging.info("common counter={}".format(counter))
    logging.info("counter_empty={}".format(counter_empty))
    logging.info("right_rst_counter={}".format(right_rst_counter))
    logging.info("rst_dict:\n{}\nRst_dict by species folder:".format(rst_dict))
    for key, value in rst_dict.items():
        logging.info("{}:\n{}".format(key, value))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with .fas files to sort', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder)
    except:
        logging.exception("Unexpected error")
