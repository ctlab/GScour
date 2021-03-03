#!/usr/bin/env python
import scandir
import argparse
import sys
import logging
import os
import logging


log_file = os.path.join(os.getcwd(), 'swamp_log.log')
print log_file
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=log_file)
BROKEN_FILES = list()


def parse_dir(infolder):
    for personal_folder in scandir.scandir(infolder):
        full_path =  os.path.join(infolder, personal_folder.name)
        if os.path.isdir(full_path):
            for infile in scandir.scandir(full_path):
                if infile.name.split('.')[-1] == 'phy':
                    yield full_path, infile


def run_swamp(personal_folder, executable_path, infile, branchcodes, threshold, windowsize):
    global BROKEN_FILES
    try:
        launch_swamp = 'python2 {}' \
                       ' -i {} -b {} -t {} -w {} >> {}'.format(executable_path, personal_folder, branchcodes,
                                                               threshold, windowsize, log_file)
        os.system(launch_swamp)
    except ValueError:
        logging.exception("infile {}: {}".format(infile, sys.exc_info()))
        BROKEN_FILES.append(infile.name.split('.')[0])
        
        



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', help='Path to the SWAMP executable', nargs='?')
    parser.add_argument('-i', help='INFOLDER', nargs='?')
    parser.add_argument('-b', help='BRANCHNAMESFILE', nargs='?')
    parser.add_argument('-t', help='THRESHOLD', nargs='?', default="")
    parser.add_argument('-w', help='WINDOWSIZE', nargs='?')
    args = parser.parse_args()
    try:
        for folder, infile in parse_dir(args.i):
            print "pass dir", folder
            run_swamp(folder, args.e, infile, args.b, args.t, args.w)
    except:
        logging.exception("Unexpected error")
    logging.info("BROKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    logging.info("The work has been completed")
