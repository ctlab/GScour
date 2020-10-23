#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import sys

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
ALIGNED_FILES = list()
EXCEPTION_NUMBER = 0

def parse_dir(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fna':
            yield os.path.join(infolder, infile)


def launch_prank(infile, outfolder, tree):
    outfile_path_without_extension = os.path.join(outfolder, re.search(r'\/(\d+)\.', infile).group(1))
    file_number = re.search(r'\/(\d+)\.', infile).group(1)
    final_file_path = '{}.best.fas'.format(outfile_path_without_extension)
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile_path_without_extension), 'log'))
    try:
        global ALIGNED_FILES
        if os.path.isfile(final_file_path):
            raise Exception('final_file_path {} already exists'.format(final_file_path))
        if tree:
            launch = 'prank -d={0} -o={1} -t={2} > {3}'.format(infile, outfile_path_without_extension, tree, log_file)
        else:
            launch = 'prank -d={0} -o={1} -showtree > {2}'.format(infile, outfile_path_without_extension, log_file)
        if not os.system(launch):
            logging.info("prank completed task for file {}".format(file_number))
            if file_number not in ALIGNED_FILES:
                ALIGNED_FILES.append(file_number)
    except:
        global EXCEPTION_NUMBER
        logging.exception("sys.exc_info() {0}, outfile number {1}".format(sys.exc_info(),
                                                                          outfile_path_without_extension))
        EXCEPTION_NUMBER += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--outfolder', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--tree', help='Path to the tree', nargs='?', default="")
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(threads)
        inputs = list(parse_dir(args.infolder))
        len_inputs = len(inputs)
        outfolder = args.outfolder
        if not os.path.isdir(outfolder):
            os.makedirs(outfolder)
        pool.starmap(launch_prank, zip(inputs, len_inputs * [outfolder], len_inputs * [args.tree]))
    except:
        logging.exception("Unexpected error")
        logging.info("Number of ALIGNED_FILES = {}".format(ALIGNED_FILES))
        logging.info("Number of prank exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("Number of ALIGNED_FILES = {}".format(ALIGNED_FILES))
    logging.info("Number of prank exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")
