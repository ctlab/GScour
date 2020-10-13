#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import time
import datetime

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
list_of_files_numbers = list()


def parse_dir(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fas':
            list_of_files_numbers.append(int(infile.split('.')[0]))
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
    launch = '/home/alina_grf/BIOTOOLS/Gblocks_0.91b/Gblocks {0} -t=d -b1=3 -b2=3 -b3=7 -b4=3 -b5=h ' \
             '-p=Yes'.format(infile)
    os.system(launch)


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
        lack = list(set(list(range(1, 3561))) - set(list_of_files_numbers))
        logging.info("lack list: {}".format(lack))
        pool.map(launch_gblocks, inputs)
    except:
        logging.exception("Unexpected error")
    logging.info("The work has been completed")
