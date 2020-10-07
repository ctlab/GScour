#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        yield os.path.join(infolder, infile)


def launch_prank(infile, outfolder, tree):
    outfile = os.path.join(outfolder, re.search(r'\/(\d+)\.', infile).group(1))
    output_file = os.path.join("{}.{}".format(os.path.abspath(outfile), 'best.fas'))
    if os.path.isfile(output_file):
        logging.info("Output file {} is already exist, skip to the next".format(output_file))
        return
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile), 'txt'))
    if tree:
        launch = 'prank -d={0} -o={1} -t={2} > {3}'.format(infile, outfile, tree, log_file)
    else:
        launch = 'prank -d={0} -o={1} > {3}'.format(infile, outfile, tree, log_file)
    os.system(launch)


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
