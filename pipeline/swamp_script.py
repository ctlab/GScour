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
                    yield os.path.join(infolder, personal_folder, infile)


def run_swamp(infile):
    file_out_path = infile.replace('phy', 'out')
    personal_dir = os.path.split(file_out_path)[0]
    branchcodes = '/home/alina_grf/progprojects/search_for_positive_selection/branchcodes.txt'
    threshold = 10
    windowsize = 3
    try:
        launch_swamp = '/home/alina_grf/BIOTOOLS/SWAMP-master/SWAMP.py' \
                       ' -i {} -b {} -t {} -w {}'.format(personal_dir, branchcodes, threshold, windowsize)
        if os.system(launch_swamp):
            logging.info("SWAMP analysis has been done for file {}".format(infile))
    except:
        logging.exception("sys.exc_info() {0}, outfile {1}".format(sys.exc_info(), file_out_path))


def main(infolder, tree):
    for infile in parse_dir(infolder):
        run_swamp(infile, tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for paml', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder, args.tree)
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
