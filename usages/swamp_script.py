#!/usr/bin/env python
import argparse
import sys
import logging
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

"""
def parse_dir(infolder):
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    yield os.path.join(infolder, personal_folder, infile)
"""

def run_swamp():
    #file_out_path = infile.replace('phy', 'out')
    #personal_dir = os.path.split(file_out_path)[0]
    common_dir = '/for_paml/'
    branchcodes = '/home/alina_grf/progprojects/search_for_positive_selection/branchcodes.txt'
    threshold = 10
    windowsize = 3
    log_file = '/home/alina_grf/progprojects/search_for_positive_selection/swamp_log.log'
    try:
        launch_swamp = 'python2 /home/alina_grf/BIOTOOLS/SWAMP-master/SWAMP.py' \
                       ' -i {} -b {} -t {} -w {} >> {}'.format(common_dir, branchcodes, threshold, windowsize, log_file)
        os.system(launch_swamp)
    except:
        logging.exception("sys.exc_info() {0}, outfile {1}".format(sys.exc_info()))


run_swamp()