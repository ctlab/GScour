#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import re
import traceback

from Bio import SeqIO

"""
This script sorts fasta files from one folder to child folders with unique names 
corresponding to set of species in fasta file. For example:
$ cd parent_dir
$ ls
1.fasta 2.fasta 3.fasta
$ less 1.fasta      $ less 2.fasta      $ less 3.fasta 
>1                  >2                  >3
ATG....             ATG...              ATG...
>2                  >3                  >2    
ATG....             ATG...              ATG...

result of the script:
$ cd parent_dir
$ ls */
12/:            23/:
1.fasta         2.fasta
                3.fasta
                
So, there will be a folder for every combination of species (for every group)
"""

LOG_FILE = "sort_by_groups.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def parse_fasta_dir(in_dir):
    """.fas-gb is format out of Gblocks
    this extension can be adjusted for needs"""
    for infile in os.listdir(in_dir):
        if infile.split('.')[-1] == 'fas-gb':
            yield os.path.join(in_dir, infile)


def form_dict_of_groups(in_dir):
    groups_dict = dict()
    for fasta_file in parse_fasta_dir(in_dir):
        list_names = list()
        file_number = re.search(r'(\d+)\.', fasta_file).group(1)
        for record in SeqIO.parse(fasta_file, "fasta"):
            try:
                list_names.append(int(record.name))
            except BaseException as err:
                logging.exception("File {}: {}, \ntraceback: P{}".format(fasta_file, err.args,
                                                                         traceback.print_tb(err.__traceback__)))
                continue
        list_names.sort()
        sorted_name = ""
        for item in list_names:
            sorted_name += str(item)
        # name_set = frozenset(sorted_name)
        if not groups_dict.get(sorted_name):
            groups_dict[sorted_name] = list()
        groups_dict[sorted_name].append(file_number)
    return groups_dict


def replace_files_by_groups(folder_in, groups_dict):
    """replace files to new group folder
    the extension of fasta file can be adjusted for needs"""
    for name_set, list_of_files in groups_dict.items():
        folder_name = "".join(name_set)
        full_folder_name = os.path.join(folder_in, folder_name)
        os.makedirs(full_folder_name)
        for file_number in list_of_files:
            os.replace(os.path.join(folder_in, file_number + ".best.fas-gb"), os.path.join(full_folder_name, file_number
                                                                                           + ".best.fas-gb"))


def main(folder_in):
    logging.info("folder to work with: {}".format(folder_in))
    groups_dict = form_dict_of_groups(folder_in)
    replace_files_by_groups(folder_in, groups_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with .fas files to sort', nargs='?')
    args = parser.parse_args()
    try:
        main(args.i)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
    logging.info("The work has been completed")
