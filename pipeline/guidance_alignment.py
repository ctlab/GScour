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
        if infile.split('.')[-1] == 'fas':
            yield os.path.join(infolder, infile)


def launch_guidance(infile, outfolder, threads):
    personal_outdir = os.path.join(outfolder, re.search(r'\/(\d+)\.', infile).group(1))
    if not os.path.isdir(personal_outdir):
        os.makedirs(personal_outdir)
    # launch = 'perl guidance.pl --seqFile {0} --msaProgram PRANK --seqType nuc --msaFile {0} --proc_num {1}' \
    # ' --outDir {2}'.format(infile, threads, outdir)
    launch = 'perl /home/alina_grf/BIOTOOLS/guidance.v2.02/www/Guidance/guidance.pl --seqFile {0} --msaProgram PRANK ' \
             '--seqType nuc ' \
             '--proc_num {1} --outDir {2} --bootstraps 20'.format(infile, threads, personal_outdir)
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
        inputs = list(parse_dir(args.infolder))
        len_inputs = len(inputs)
        pool.starmap(launch_guidance, zip(inputs, len_inputs * [args.outfolder], len_inputs * [args.threads]))
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
