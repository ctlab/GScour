#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
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


def parse_dir(in_dir):
    for infile in os.listdir(in_dir):
        if infile.split('.')[-1] == 'fas':
            yield os.path.join(in_dir, infile)


def launch_gblocks(infile, exec_path):
    global counter_file
    file_number = re.search(r'\/(\d+)\.', infile).group(1)
    launch = '{} {} -t=c -b1=5 -b2=5 -b3=7 -b4=2 -b5=h ' \
             '-p=Yes >> {}'.format(exec_path, infile, LOG_FILE)  # TODO: > LOG_FILE: to do multiprocessing
    os.system(launch)
    logging.info("Gblocks processed file {} with params {}".format(file_number, '-t=c -b1=3 -b2=4 -b3=7 -b4=6 -b5=h'))
    with counter_file.get_lock():
        counter_file.value += 1
        logging.info("Counter (processed files) = {}".format(counter_file.value))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with input files for Gblocks'
                                    'FASTA formats are accepted', nargs='?')
    parser.add_argument('--exec', help='Path to the Gblocks executable', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        counter_file = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter_file,))
        in_dir = args.i
        executable_path = args.exec
        inputs = list(parse_dir(in_dir))
        len_inputs = len(inputs)
        logging.info("Path to the folder with input files for Gblocks: {}\nExecutable path: {}".
                     format(in_dir, executable_path))
        i = pool.starmap_async(launch_gblocks, zip(inputs, len_inputs * [executable_path]))
        i.wait()
        i.get()
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        logging.info("Number of processed files = {}".format(counter_file.value))

    logging.info("Number of processed files = {}".format(counter_file.value))
    logging.info("The work has been completed")
