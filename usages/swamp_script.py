#!/usr/bin/env python
import argparse
import sys
import logging
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def parse_dir(infolder):
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    yield os.path.join(infolder, personal_folder)


def run_swamp(personal_folder, branchcodes, threshold, windowsize):
    log_file = 'swamp_log.log'
    try:
        launch_swamp = 'python2 /home/alina_grf/BIOTOOLS/SWAMP-master/SWAMP.py' \
                       ' -i {} -b {} -t {} -w {} >> {}'.format(personal_folder, branchcodes,
                                                               threshold, windowsize, log_file)
        os.system(launch_swamp)
    except ValueError:
        logging.exception("sys.exc_info() {0}".format(sys.exc_info()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('-b', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('-t', help='Path to the tree', nargs='?', default="")
    parser.add_argument('-w', help='Number of threads', nargs='?')
    args = parser.parse_args()
    try:
        for folder in parse_dir(args.infodler):
            run_swamp(folder, args.branchcodes, args.threshold, args.windowsize)
    except:
        logging.exception("Unexpected error")
    logging.info("The work has been completed")