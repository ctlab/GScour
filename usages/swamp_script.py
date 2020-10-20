#!/usr/bin/env python
import argparse
import re
import sys
import logging
import os
import logging

LOG_FILE = 'swamp_log.log'
BROCKEN_FILES = list()
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def parse_dir(infolder):
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    yield os.path.join(infolder, personal_folder)


def run_swamp(personal_folder, branchcodes, threshold, windowsize):
    try:
        launch_swamp = 'python2 /home/alina_grf/BIOTOOLS/SWAMP-master/SWAMP.py' \
                       ' -i {} -b {} -t {} -w {} >> {}'.format(personal_folder, branchcodes,
                                                               threshold, windowsize, LOG_FILE)
        os.system(launch_swamp)
    except ValueError:
        logging.exception("sys.exc_info() {0}".format(sys.exc_info()))
        file_number = (re.search(r"/(\d+)/*", personal_folder)).group(1)
        if file_number not in BROCKEN_FILES:
            BROCKEN_FILES.append(file_number)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('-b', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('-t', help='Path to the tree', nargs='?', default="")
    parser.add_argument('-w', help='Number of threads', nargs='?')
    args = parser.parse_args()
    try:
        for folder in parse_dir(args.i):
            run_swamp(folder, args.b, args.t, args.w)
        if BROCKEN_FILES:
            logging.warning("BROCKEN_FILES: {}".format(BROCKEN_FILES))
    except:
        logging.exception("Unexpected error")
    logging.info("The work has been completed")