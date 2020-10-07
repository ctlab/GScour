#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import multiprocessing
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        yield os.path.join(infolder, infile)


def rename_seq(infile, outfolder_fas, outfolder_txt):
    if not os.path.isdir(outfolder_fas):
        os.makedirs(outfolder_fas)
    if not os.path.isdir(outfolder_txt):
        os.makedirs(outfolder_txt)

    with open(infile, "r+") as handle:
        for index, record in enumerate(SeqIO.parse(handle, "fasta")):
            file_identifier = re.search(r'\/(\d+)\.', infile).group(1)
            with open(os.path.join(outfolder_txt, file_identifier + ".txt"), "a") as ofile:
                logging.info("writing file with gene id {}".format(ofile.name))
                ofile.write("{}\n".format(record.id))
            seq = SeqRecord(record.seq, id=str(index+1), description="")
            with open(os.path.join(outfolder_fas, file_identifier + ".fna"), "a") as ofile:
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
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")