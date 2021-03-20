#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import traceback

PROCESSED_FILES = list()
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
    global PROCESSED_FILES
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
     [please provide full (and not relative) path]
    --bootstraps: Number of bootstrap iterations (only for GUIDANCE). Default=100
    --prank: path to prank executable. Default=prank
    """
    launch = 'perl {} --seqFile {} --msaProgram PRANK ' \
             '--seqType nuc ' \
             '--proc_num {} --outDir {} --bootstraps 20'.format(executable_path, infile, number_of_threads,
                                                                personal_dir_out)

    os.system(launch)
    logging.info("Guidance completed task for file {}".format(file_number))
    if file_number not in PROCESSED_FILES:
        PROCESSED_FILES.append(file_number)  # TODO: to shared variables
        with counter.get_lock():
            counter.value += 1
            logging.info("Counter (PROCESSED_FILES) = {}\nList of PROCESSED_FILES: {}".
                         format(counter.value, PROCESSED_FILES))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with input files (.fna) for guidance', nargs='?')
    parser.add_argument('--o', help='Path to the folder with output files of guidance', nargs='?')
    parser.add_argument('--exec', help='Path to the guidance executable "guidance.pl"', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    in_dir = args.i
    out_dir = args.o
    logging.info("Path to the folder with input files for guidance: {}\n"
                 "Path to the folder with output files of guidance: {}\n".format(in_dir, out_dir))
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        inputs = list(parse_dir(in_dir))
        len_inputs = len(inputs)
        i = pool.starmap_async(launch_guidance, zip(inputs, len_inputs * [out_dir], len_inputs * [threads],
                                                    len_inputs * [args.exec]))
        i.wait()
        i.get()
    except BaseException as e:
        logging.exception("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))

    logging.info("Number of ALIGNED_FILES = {}".format(counter.value))
    logging.info("Number of exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")