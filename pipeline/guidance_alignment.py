#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import sys

ALIGNED_FILES = list()
EXCEPTION_NUMBER = 0
counter = None
LOG_FILE = "guidance_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(args_counter):
    """store the counter for later use"""
    global counter
    counter = args_counter


def parse_dir(folder_in):
    for infile in os.listdir(folder_in):
        if infile.split('.')[-1] == 'fna':
            yield os.path.join(folder_in, infile)


def launch_guidance(infile, folder_out, number_of_threads, executable_path):
    global counter
    global ALIGNED_FILES
    global EXCEPTION_NUMBER
    file_number = re.search(r'\/(\d+)\.', infile).group(1)
    personal_dir_out = os.path.join(folder_out, file_number)
    if not os.path.isdir(personal_dir_out):
        os.makedirs(personal_dir_out)
    """
    Required parameters:
    --seqFile: Input sequence file in FASTA format
    --msaProgram: Which MSA program to use
    --seqType: Type of sequences for alignment (amino acids, nucleotides, or codons)
    --outDir: Output directory that will be created automatically and hold all output files
     [please provid full (and not relative) path]
    --bootstraps: Number of bootstrap iterations (only for GUIDQANCE). Defaut=100
    --prank: path to prank executable. Default=prank
    """
    launch = 'perl {} --seqFile {} --msaProgram PRANK ' \
             '--seqType nuc ' \
             '--proc_num {} --outDir {} --bootstraps 20'.format(executable_path, infile, number_of_threads,
                                                                personal_dir_out)
    try:
        if not os.system(launch):
            logging.info("guidance completed task for file {}".format(file_number))
            if file_number not in ALIGNED_FILES:
                ALIGNED_FILES.append(file_number)
                with counter.get_lock():
                    counter.value += 1
                    logging.info("Counter (ALIGNED_FILES) = {}\nList of ALIGNED_FILES: {}".
                                 format(counter.value, ALIGNED_FILES))
    except:  # TODO: don't catch need to be fixed
        logging.exception("{}, file_number {}".format(sys.exc_info(), file_number))
        EXCEPTION_NUMBER += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--outfolder', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--exec', help='Path to the guidance executable "guidance.pl"', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    infolder = args.infolder
    outfolder = args.outfolder
    logging.info("Path to the folder with input files for guidance: {}\n"
                 "Path to the folder with output files of guidance: {}\n".format(infolder, outfolder))
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        inputs = list(parse_dir(infolder))
        len_inputs = len(inputs)
        i = pool.starmap_async(launch_guidance, zip(inputs, len_inputs * [outfolder], len_inputs * [threads], len_inputs * [
            args.exec]))
        i.wait()
        i.get()
    except:
        logging.exception("Unexpected error")

    logging.info("Number of ALIGNED_FILES = {}".format(counter.value))
    logging.info("Number of exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")
