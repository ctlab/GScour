#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import math
import multiprocessing
import os
import logging
import re
import traceback

import pandas as pd

counter_file = None
CORRECT_FILES = list()
ERROR_FILES = list()
CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()
LOG_FILE = "gblocks_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(args):
    """store the counters for later use"""
    global counter_file
    counter_file = args


def parse_dir(input_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fas':
                    yield input_dir, species_folder.name, infile.name


def gather_correct(correct):
    for tup in correct:
        for correct_item in tup[0]:
            CORRECT_FILES_TO_WRITE.add(correct_item)
    logging.info("CORRECT_FILES_TO_WRITE {} {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))


def gather_errors(exception):
    for tup in exception:
        for error_item in tup[1]:
            ERROR_FILES_TO_WRITE.add(error_item)
    logging.info("ERROR_FILES_TO_WRITE {} {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))


def launch_gblocks(input_tuple, auto_flag, exec_path, child_logger):
    global counter_file
    global CORRECT_FILES
    global ERROR_FILES
    try:
        input_dir, species_folder, infile_name = input_tuple
        file_gene_name = infile_name.split('.')[0]
        infile_path = os.path.join(input_dir, species_folder, infile_name)
        if auto_flag == 'n':
            params_string = '-t=c -b1=3 -b2=4 -b3=8 -b4=9 -b5=n -p=y'
            logging.info("auto flag = 'n', custom parameter string is\n{}".format(params_string))
        else:
            if len(species_folder) <= 9:
                number_of_species = len(species_folder)
            else:
                number_of_species = (len(species_folder) - 9) / 2 + 9
            b1 = math.ceil(number_of_species / 2 + 1)
            b2 = math.ceil(number_of_species * 0.85)
            params_string = '-t=c -b1={} -b2={} -b3=8 -b4=9 -b5=n -p=y'.format(b1, b2)
            logging.info("auto flag = 'y', auto calculating parameter string is\n{}".format(params_string))

        launch = '{} {} {} >> {}'.format(exec_path, infile_path, params_string, LOG_FILE)
        os.system(launch)
        child_logger.info("Gblocks processed file {} with params {}".format(infile_name, params_string))
        CORRECT_FILES.append(file_gene_name)
        with counter_file.get_lock():
            counter_file.value += 1
            child_logger.info("Counter (processed files) = {}".format(counter_file.value))
    except BaseException as err:
        ERROR_FILES.append(file_gene_name)
        logging.exception("Infile {}, - Unexpected error: {}".format(infile_name, err))
    return CORRECT_FILES, ERROR_FILES


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for Gblocks'
                                    'FASTA formats are accepted', nargs='?', required=True)
    parser.add_argument('--exec', help='Path to the Gblocks executable', nargs='?', required=True)
    parser.add_argument('--auto', help='\'y\' or \'n\': automatic selection of basic parameters according to group size',
                        nargs='?', required=True)
    parser.add_argument('--threads', help='Number of threads', nargs='?', required=True)
    args = parser.parse_args()
    threads = int(args.threads)
    try:
        counter_file = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = logging.getLogger(__name__)
        logger.addHandler(logging.FileHandler(LOG_FILE))
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter_file,))
        in_dir = args.i
        executable_path = args.exec
        input_tuples = list(parse_dir(in_dir))
        len_inputs = len(input_tuples)
        logger.info("Path to the folder with input files for Gblocks: {}\nExecutable path: {}".
                    format(in_dir, executable_path))
        i = pool.starmap_async(launch_gblocks, zip(input_tuples, len_inputs * [args.auto],
                                                   len_inputs * [executable_path], len_inputs * [logger]),
                               callback=gather_correct,
                               error_callback=gather_errors)
        i.wait()
        i.get()
    except BaseException as e:
        logger.exception("Unexpected error: {}".format(e))
        logger.info("Number of processed files = {}".format(counter_file.value))
        raise e
    # logger.info("Number of processed files = {}".format(counter_file.value))
    logging.info("CORRECT_FILES_TO_WRITE {} {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))
    logging.info("ERROR_FILES_TO_WRITE {} {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))
    resulting_file = os.path.join(args.i, 'gblocks_summary.xlsx')
    logging.info("res file {}".format(resulting_file))
    writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
    df_corr = pd.DataFrame({'Gene name': list(CORRECT_FILES_TO_WRITE)})
    df_corr.to_excel(writer, sheet_name='correct files', index=False)
    df_err = pd.DataFrame({'Gene name': list(ERROR_FILES_TO_WRITE)})
    df_err.to_excel(writer, sheet_name='exception files', index=False)
    writer.save()
    logger.info("The work has been completed")
