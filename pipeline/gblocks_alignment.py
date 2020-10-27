#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import sys


logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
ALIGNED_FILES = list()
EXCEPTION_NUMBER = 0


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fas':
            yield os.path.join(infolder, infile)


def launch_gblocks(infile):
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

    global ALIGNED_FILES
    file_number = re.search(r'\/(\d+)\.', infile).group(1) #/home/alina_grf/BIOTOOLS/Gblocks_0.91b/Gblocks
    launch = '/home/alina_grf/BIOTOOLS/Gblocks_0.91b/Gblocks {0} -t=d -b1=3 -b2=0 -b3=7 -b4=6 -b5=h ' \
             '-p=Yes'.format(infile)
    if not os.system(launch):
        logging.info("gblocks completed task for file {}".format(file_number))
        if file_number not in ALIGNED_FILES:
            ALIGNED_FILES.append(file_number)
            print(ALIGNED_FILES)
    else:
        global EXCEPTION_NUMBER
        EXCEPTION_NUMBER += 1
        print(EXCEPTION_NUMBER)
        logging.exception("gblocks error for file number {}".format(file_number))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(threads)
        inputs = list(parse_dir(args.infolder))
        pool.map(launch_gblocks, inputs)
    except:
        logging.info("Number of ALIGNED_FILES = {}".format(ALIGNED_FILES))
        logging.info("Number of gblocks exceptions = {}".format(EXCEPTION_NUMBER))
        logging.exception("Unexpected error")
    global ALIGNED_FILES
    global EXCEPTION_NUMBER
    print("al", ALIGNED_FILES)
    print("excep", EXCEPTION_NUMBER)
    logging.info("Number of ALIGNED_FILES = {}".format(ALIGNED_FILES))
    logging.info("Number of gblocks exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")
