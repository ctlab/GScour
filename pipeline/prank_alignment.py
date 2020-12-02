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
LOG_FILE = "prank_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(counter_args):
    """store the counter for later use"""
    global counter
    counter = counter_args


def parse_dir(folder_in):
    for infile in os.listdir(folder_in):
        if infile.split('.')[-1] == 'fna':
            yield os.path.join(folder_in, infile)


def get_final_file_path(outfile_path_without_extension, format_out):
    if format_out in ["phylipi", "phylips", "paml"]:
        final_file_path = '{}.best.phy'.format(outfile_path_without_extension)
    elif format_out in ["fasta", ""]:
        final_file_path = '{}.best.fas'.format(outfile_path_without_extension)
    else:
        final_file_path = '{}.best.nex'.format(outfile_path_without_extension)
    return final_file_path


def get_launch_command(infile, final_file_path, outfile_path_without_extension, tree, format_out):
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile_path_without_extension), 'log'))
    if os.path.isfile(final_file_path):
        raise Exception('final_file_path {} already exists'.format(final_file_path))
    if tree and format_out:
        launch = 'prank -d={0} -o={1} -t={2} -codon -f={3} > {4}'.format(infile, outfile_path_without_extension,
                                                                         tree, format_out, log_file)
    elif tree:
        launch = 'prank -d={0} -o={1} -t={2} -codon > {3}'.format(infile, outfile_path_without_extension,
                                                                  tree, log_file)
    elif format_out and not tree:
        launch = 'prank -d={0} -o={1} -showtree -codon -f=paml > {2}'.format(infile,
                                                                             outfile_path_without_extension,
                                                                             log_file)
    else:
        launch = 'prank -d={0} -o={1} -showtree -codon  > {2}'.format(infile,
                                                                      outfile_path_without_extension,
                                                                      log_file)
    return launch


def launch_prank(infile, folder_out, tree, format_out):
    global counter
    outfile_path_without_extension = os.path.join(folder_out, re.search(r'\/(\d+)\.', infile).group(1))
    file_number = re.search(r'\/(\d+)\.', infile).group(1)
    final_file_path = get_final_file_path(outfile_path_without_extension, format_out)
    try:
        global ALIGNED_FILES
        launch_command = get_launch_command(infile, final_file_path, outfile_path_without_extension, tree, format_out)
        if not os.system(launch_command):
            logging.info("prank completed task for file {}".format(file_number))
            if file_number not in ALIGNED_FILES:
                ALIGNED_FILES.append(file_number)
                with counter.get_lock():
                    counter.value += 1
                    logging.info("Counter (ALIGNED_FILES) = {}".format(counter.value))

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
    parser.add_argument('--f', help='Output format: ["fasta" (default), "phylipi", "phylips", "paml", "nexus"]',
                        nargs='?', default="")
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        infolder = args.infolder
        inputs = list(parse_dir(args.infolder))
        len_inputs = len(inputs)
        outfolder = args.outfolder
        logging.info("Path to the folder with input files for prank: {}\n"
                     "Path to the folder with output files of prank: {}".format(infolder, outfolder))
        if not os.path.isdir(outfolder):
            os.makedirs(outfolder)
        output_format = args.f
        if output_format not in ["fasta", "phylipi", "phylips", "paml", "nexus", ""]:
            raise SyntaxError("Not valid output format, check option --f, -h for help")
        i = pool.starmap_async(launch_prank, zip(inputs, len_inputs * [outfolder], len_inputs * [args.tree],
                                                 len_inputs * [output_format]))
        i.wait()
        i.get()
        # pool.starmap(launch_prank, zip(inputs, len_inputs * [outfolder], len_inputs * [args.tree]))
    except:
        logging.exception("Unexpected error")
        logging.info("Number of ALIGNED_FILES = {}".format(counter.value))
        logging.info("Number of prank exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("Number of ALIGNED_FILES = {}".format(counter.value))
    logging.info("Number of prank exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")
