#!/usr/bin/env python
import argparse
import re
import sys
import logging
import os
import logging

LOG_FILE = 'swamp_log.log'
BROKEN_FILES = list()
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def parse_dir(infolder):
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    yield os.path.join(infolder, personal_folder)


def run_swamp(personal_folder, exec_path, branchcodes, threshold, windowsize):
    global BROKEN_FILES
    try:
        launch_swamp = '{}' \
                       ' -i {} -b {} -t {} -w {} >> {}'.format(exec_path, personal_folder, branchcodes,
                                                               threshold, windowsize, LOG_FILE)
        if os.system(launch_swamp):
            raise ValueError
    except ValueError:
        logging.exception("sys.exc_info() {0}".format(sys.exc_info()))
        file_number = (re.search(r"/(\d+)/*", personal_folder)).group(1)
        if file_number not in BROKEN_FILES:
            BROKEN_FILES.append(file_number)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-exec', help='Path to the SWAMP executable', nargs='?')
    parser.add_argument('-i', help='INFOLDER', nargs='?')
    parser.add_argument('-b', help='BRANCHNAMESFILE', nargs='?')
    parser.add_argument('-t', help='THRESHOLD', nargs='?', default="")
    parser.add_argument('-w', help='WINDOWSIZE', nargs='?')
    args = parser.parse_args()
    try:
        for folder in parse_dir(args.i):
            run_swamp(folder, args.exec, args.b, args.t, args.w)
    except:
        logging.exception("Unexpected error")
    logging.info("BROKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    logging.info("The work has been completed")