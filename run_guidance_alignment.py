#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        yield os.path.join(infolder, infile)


def launch_guidance(infile, outfolder):
    outdir = os.path.join(outfolder, os.path.splitext(os.path.splitext(infile)[0])[0])
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # launch = 'perl guidance.pl --seqFile {0} --msaProgram PRANK --seqType nuc --msaFile {0} --proc_num {1}' \
    # ' --outDir {2}'.format(infile, threads, outdir)
    launch = 'perl guidance.pl --seqFile {0} --seqType nuc --proc_num {1} --outDir {2}'.format(infile, threads, outdir)
    os.system(launch)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--outfolder', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(threads)
        pool.starmap(launch_guidance, zip(parse_dir(args.infolder), args.outfolder, threads))
    except:
        logging.exception("Unexpected error")
