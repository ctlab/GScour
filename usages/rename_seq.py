#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import multiprocessing
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
""" Rename sequences in fasta file: change name to number"""


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        yield os.path.join(infolder, infile)


def rename_seq(infile, out_folder_fas, out_folder_log):
    if not os.path.isdir(out_folder_fas):
        os.makedirs(out_folder_fas)
    if not os.path.isdir(out_folder_log):
        os.makedirs(out_folder_log)

    with open(infile, "r+") as handle:
        for index, record in enumerate(SeqIO.parse(handle, "fasta")):
            file_identifier = re.search(r'\/(\d+)\.', infile).group(1)
            with open(os.path.join(out_folder_log, file_identifier + ".log"), "a") as ofile:
                logging.info("writing file with gene id {}".format(ofile.name))
                ofile.write("{}\n".format(record.id))
            seq = SeqRecord(record.seq, id=str(index+1), description="")
            with open(os.path.join(out_folder_fas, file_identifier + ".fna"), "a") as ofile:
                SeqIO.write(seq, ofile, "fasta")
                logging.info("writing new fasta renaming seq file {}".format(ofile.name))
    os.remove(infile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--outfolder_fas', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--outfolder_txt', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    try:
        threads = int(args.threads)
        inputs = list(parse_dir(args.infolder))
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(threads)
        pool.starmap(rename_seq, zip(inputs, len(inputs) * [args.outfolder_fas], len(inputs) * [args.outfolder_txt]))
    except BaseException as err:
        logging.info("Unexpected error: {}".format(err))

    logging.info("The work has been completed")