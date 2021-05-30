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
LOG_FILE = "prank_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(counter_args):
    """store the counter for later use"""
    global counter
    counter = counter_args


def parse_dir(input_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fna':
                    yield input_dir, species_folder.name, infile.name


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
    if tree and format_out:  # -translate (standard code); codon alignment with the option -codon (in -codon case
        # be careful about not multiple of three sequences)
        launch = 'prank -d={0} -o={1} -t={2} -translate -f={3} > {4}'.format(infile, outfile_path_without_extension,
                                                                             tree, format_out, log_file)
    elif tree:
        launch = 'prank -d={0} -o={1} -t={2} -translate > {3}'.format(infile, outfile_path_without_extension,
                                                                      tree, log_file)
    elif format_out and not tree:
        launch = 'prank -d={0} -o={1} -showtree -translate -f={2} > {3}'.format(infile,
                                                                                outfile_path_without_extension,
                                                                                format_out,
                                                                                log_file)
    else:
        launch = 'prank -d={0} -o={1} -showtree -translate > {2}'.format(infile,
                                                                         outfile_path_without_extension,
                                                                         log_file)
    return launch


def launch_prank(input_tuple, folder_out, tree, format_out):
    global counter
    input_dir, species_folder, infile = input_tuple
    infile_path = os.path.join(input_dir, species_folder, infile)
    outfile_path_without_extension = os.path.join(folder_out, re.search(r'(\d+)\.', infile).group(1))
    file_number = re.search(r'\/(\d+)\.', infile).group(1)
    final_file_path = get_final_file_path(outfile_path_without_extension, format_out)
    try:
        global PROCESSED_FILES
        launch_command = get_launch_command(infile_path, final_file_path, outfile_path_without_extension, tree, format_out)
        if not os.system(launch_command):
            # logging.info("prank completed task for file {}".format(file_number))
            if file_number not in PROCESSED_FILES:
                PROCESSED_FILES.append(file_number)
                with counter.get_lock():
                    counter.value += 1  # TODO: wrong number
                    # logging.info("Counter (ALIGNED_FILES) = {}".format(counter.value)) # TODO: multiprocessing
    except BaseException as err:
        global EXCEPTION_NUMBER
        logging.exception("Infile {}, - Unexpected error: {}".format(infile, err))
        EXCEPTION_NUMBER += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files (.fna) for prank',
                        nargs='?')
    parser.add_argument('--o', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--tree', help='Path to the tree, exclude if there is no tree', nargs='?', default="")
    parser.add_argument('--f', help='Output format: ["fasta" (default, exclude option --f if left by default),'
                                    '"phylipi", "phylips", "paml", "nexus"]', nargs='?', default="")
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        in_dir = args.i
        input_tuples = list(parse_dir(in_dir))
        len_inputs = len(input_tuples)
        out_dir = args.o
        logging.info("Path to the folder with input files for prank: {}\n"
                     "Path to the folder with output files of prank: {}".format(in_dir, out_dir))
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        output_format = args.f
        if output_format not in ["fasta", "phylipi", "phylips", "paml", "nexus", ""]:
            raise SyntaxError("Not valid output format, check option --f, -h for help")
        i = pool.starmap_async(launch_prank, zip(input_tuples, len_inputs * [out_dir], len_inputs * [args.tree],
                                                 len_inputs * [output_format]))
        i.wait()
        i.get()
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    logging.info("The work has been completed")
