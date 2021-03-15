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


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fas':
            yield os.path.join(infolder, infile)


def launch_gblocks(infile, exec_path):
    """outdated inspection
    if os.path.splitext(os.path.splitext(infile)[1])[0] == '.fas':
        delta = datetime.datetime.now() - datetime.datetime.strptime(time.ctime(os.path.getmtime(infile)),
                                                                         "%a %b %d %H:%M:%S %Y")
        if delta.seconds > 7000:
            logging.info("The last touch time of file {} is more then 7000s".format(infile))
            launch = '/home/alina_grf/BIOTOOLS/Gblocks_0.91b/Gblocks {0} -t=d -b1=3 -b2=3 -b3=7 -b4=3 -b5=h -b6=Yes ' \
                     '-p=Yes'.format(infile)
            os.system(launch)
            """

    global counter_file
    file_number = re.search(r'\/(\d+)\.', infile).group(1)
    launch = '{} {} -t=c -b1=5 -b2=5 -b3=7 -b4=6 -b5=h ' \
             '-p=Yes'.format(exec_path, infile)              # TODO: > LOG_FILE: to do multiprocessing
    os.system(launch)
    logging.info("Gblocks processed file {} with params {}".format(file_number, '-t=c -b1=3 -b2=4 -b3=7 -b4=6 -b5=h'))
    with counter_file.get_lock():
        counter_file.value += 1
        logging.info("Counter (processed files) = {}".format(counter_file.value))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exec', help='Path to the Gblocks executable', nargs='?')
    parser.add_argument('--infolder', help='Path to the folder with input files for gblocks'
                                           'At the moment only the NBRF/PIR and FASTA formats are accepted', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        counter_file = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter_file,))
        infolder = args.infolder
        executable_path = args.exec
        inputs = list(parse_dir(infolder))
        len_inputs = len(inputs)
        logging.info("Path to the folder with input files for Gblocks: {}\nExecutable path: {}".
                     format(infolder, executable_path))
        i = pool.starmap_async(launch_gblocks, zip(inputs, len_inputs * [executable_path]))
        i.wait()
        i.get()
    except BaseException as e:
        logging.info("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))
        logging.info("Number of processed files = {}".format(counter_file.value))

    logging.info("Number of processed files = {}".format(counter_file.value))
    logging.info("The work has been completed")

