#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import math
import multiprocessing
import os
import logging
import re
import traceback

counter_file = None
LOG_FILE = "gblocks_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(args):
    """store the counters for later use"""
    global counter_file
    counter_file = args


def parse_dir(input_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fas':
                    yield input_dir, species_folder.name, infile.name


def launch_gblocks(input_tuple, auto_flag, exec_path, child_logger):
    global counter_file
    input_dir, species_folder, infile = input_tuple
    infile_path = os.path.join(input_dir, species_folder, infile)
    if auto_flag == 'n':
        params_string = '-t=c -b1=3 -b2=4 -b3=8 -b4=10 -b5=n -p=y'
    else:
        if len(species_folder) <= 9:
            number_of_species = len(species_folder)
        else:
            number_of_species = (len(species_folder) - 9) / 2 + 9
        b1 = math.ceil(number_of_species / 2 + 1)
        b2 = math.ceil(number_of_species * 0.85)
        params_string = '-t=c -b1={} -b2={} -b3=8 -b4=9 -b5=n -p=y'.format(b1, b2)

    launch = '{} {} {} >> {}'.format(exec_path, infile_path, params_string, LOG_FILE)  # TODO: > LOG_FILE: to do multiprocessing
    os.system(launch)
    child_logger.info("Gblocks processed file {} with params {}".format(infile, params_string))
    with counter_file.get_lock():
        counter_file.value += 1
        child_logger.info("Counter (processed files) = {}".format(counter_file.value))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for Gblocks'
                                    'FASTA formats are accepted', nargs='?')
    parser.add_argument('--exec', help='Path to the Gblocks executable', nargs='?')
    parser.add_argument('--auto', help='\'y\' or \'n\': automatic selection of basic parameters according to group size',
                        nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        counter_file = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = logging.getLogger(__name__)
        logger.addHandler(logging.FileHandler(LOG_FILE))
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter_file,))
        in_dir = args.i
        executable_path = args.exec
        input_tuples = list(parse_dir(in_dir))
        len_inputs = len(input_tuples)
        logger.info("Path to the folder with input files for Gblocks: {}\nExecutable path: {}".
                    format(in_dir, executable_path))
        i = pool.starmap_async(launch_gblocks, zip(input_tuples, len_inputs * [args.auto],
                                                   len_inputs * [executable_path], len_inputs * [logger]))
        i.wait()
        i.get()
    except BaseException as e:
        logger.exception("Unexpected error: {}".format(e))
        logger.info("Number of processed files = {}".format(counter_file.value))
        raise e
    logger.info("Number of processed files = {}".format(counter_file.value))
    logger.info("The work has been completed")
