#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import sys

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fna':
            yield os.path.join(infolder, infile)


def launch_prank(infile, outfolder, tree):
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    outfile_path_without_extension = os.path.join(outfolder, re.search(r'\/(\d+)\.', infile).group(1))
    final_file_path = '{}.best.fas'.format(outfile_path_without_extension)
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile_path_without_extension), 'log'))
    try:
        if os.path.isfile(final_file_path):
            raise Exception('final_file_path {} already exists'.format(final_file_path))
        if tree:
            launch = 'prank -d={0} -o={1} -t={2} > {3}'.format(infile, outfile_path_without_extension, tree, log_file)
        else:
            launch = 'prank -d={0} -o={1} -showtree > {2}'.format(infile, outfile_path_without_extension, log_file)
        os.system(launch)
    except:
        logging.exception("sys.exc_info() {0}, outfile number {1}".format(sys.exc_info(),
                                                                          outfile_path_without_extension))


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
        pool.starmap(launch_prank, zip(inputs, len_inputs * [args.outfolder], len_inputs * [args.tree]))
    except:
        logging.exception("Unexpected error")
    logging.info("The work has been completed")
