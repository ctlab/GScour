#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import re
import traceback

from Bio import SeqIO

"""
This script sorts FASTA files from one folder into child folders with names equal to the names of the FASTA files.
"""


def parse_dir_with_fasta(in_dir):
    """adjust for needs"""
    for infile in os.listdir(in_dir):
        if infile.split('.')[-1] == 'tree':
            yield in_dir, infile


def parse_dir_with_folders(in_dir):
    """ with folders and .order"""
    for species_dir in os.scandir(in_dir):
        if os.path.isdir(species_dir):
            # for infile in os.scandir(species_dir):
            #     if infile.name.split('.')[-1] == 'order':
            yield in_dir, species_dir.name  # , infile.name


def main(folder_in, extension, folder_out):
    logging.info("folder to work with: {}".format(folder_in))
    for in_dir, species_dir in parse_dir_with_folders(folder_in):
        for item in os.listdir(os.path.join(in_dir, species_dir)):
            print('item', item)
            try:
                if item.split('.')[-1] == extension:
                    new_path = '{}/{}/{}/'.format(folder_out, species_dir, item.split('.')[0])
                    os.makedirs(new_path)
                    os.popen('cp {} {}'.format(os.path.join(in_dir, species_dir, item), new_path))
            except:
                continue


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with files to sort', nargs='?')
    parser.add_argument('--e', help='Extension of input files', nargs='?')
    parser.add_argument('--o', help='Path to the output folder', nargs='?')
    args = parser.parse_args()
    if not os.path.isdir(args.o):
        os.makedirs(args.o)
    try:
        main(args.i, args.e, args.o)
    except BaseException as e:
        logging.info("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))
    logging.info("The work has been completed")
