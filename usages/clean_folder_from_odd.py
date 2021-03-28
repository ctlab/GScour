#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os

BROKEN_FASTA = list()
REMOVED_FASTA = list()
REMOVED_LOG = list()
NOT_FOUND_FILES = list()
LOG_FILE = "clean_folder_from_odd.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def parse_broken_dir(infolder):
    global BROKEN_FASTA
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fna':
            BROKEN_FASTA.append(infile.split('.')[0])


def get_files(infolder):
    fasta_files = list()
    log_list = list()
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fas':
            fasta_files.append(infile.split('.')[0])
        elif infile.split('.')[-1] == 'log':
            log_list.append(infile.split('.')[0])
    return fasta_files, log_list


def clean_dir_from_odd(infolder):
    global BROKEN_FASTA
    global REMOVED_FASTA
    global REMOVED_LOG
    global NOT_FOUND_FILES
    NOT_FOUND_FILES = BROKEN_FASTA
    for infile in os.listdir(infolder):
        file_name = infile.split('.')[0]
        tail = infile.split('.')[-1]
        try:
            if file_name in BROKEN_FASTA:
                os.remove(os.path.join(infolder, infile))
                NOT_FOUND_FILES.remove(file_name)
                if tail == "fas":
                    REMOVED_FASTA.append(file_name)
                elif tail == "log":
                    REMOVED_LOG.append(file_name)
        except FileNotFoundError:
            if file_name not in NOT_FOUND_FILES:
                NOT_FOUND_FILES.append(file_name)


def clean_folder_from_redundant_log(folder_in):
    fasta_list, log_list = get_files(folder_in)
    residue = [x for x in log_list if x not in fasta_list]
    logging.info("residue list of redundant log files: {}".format(residue))
    try:
        for res in residue:
            os.remove(os.path.join(folder_in, "{}.log".format(res)))
    except FileNotFoundError:
        logging.warning("not found log file number {}".format(res))
    logging.info("files from residue list removed")


def main(broken_folder, cleaning_folder):
    parse_broken_dir(broken_folder)
    clean_dir_from_odd(cleaning_folder)
    clean_folder_from_redundant_log(cleaning_folder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--broken', help='Path to the folder with broken files', nargs='?')
    parser.add_argument('--clean', help='Path to the folder for cleaning', nargs='?')
    args = parser.parse_args()
    try:
        main(args.broken, args.clean)
    except BaseException as err:
        logging.info("Unexpected error: {}".format(err))
    logging.info("BROKEN_FASTA list of length {}:\n{}".format(len(BROKEN_FASTA), BROKEN_FASTA))
    logging.info("REMOVED_FASTA list of length {}:\n{}".format(len(REMOVED_FASTA), REMOVED_FASTA))
    logging.info("REMOVED_LOG list of length {}:\n{}".format(len(REMOVED_LOG), REMOVED_LOG))
    logging.info("NOT_FOUND_FILES list of length {}:\n{}".format(len(NOT_FOUND_FILES), NOT_FOUND_FILES))