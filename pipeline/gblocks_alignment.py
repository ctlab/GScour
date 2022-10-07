#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import math
import os
import logging
import sys

import pandas as pd
import functools, multiprocessing, inspect


CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()
LOG_FILE = "gblocks_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def parse_dir(input_dir, extension):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == extension:
                    yield input_dir, species_folder.name, infile.name


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except:
            frames = inspect.trace()
            argvalues = inspect.getargvalues(frames[0][0])
            file_name = argvalues.locals['args'][0][2]
            gene_name = file_name.split('.')[0]
            logging.exception("gene {}".format(gene_name))
            logging.exception("{} exception in file {}/{}".format(sys.exc_info()[0],
                                                                  argvalues.locals['args'][0][1],
                                                                  gene_name))
    return wrapped_func


@trace_unhandled_exceptions
def launch_gblocks(input_tuple, auto_flag, exec_path, child_logger):
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
    child_logger.info("Launching Gblocks for file {} with params {}".format(infile_name, params_string))
    os.system(launch)
    child_logger.info("OK for file {}".format(file_gene_name))


def get_correct_and_errors_from_log_file():
    global ERROR_FILES_TO_WRITE
    global CORRECT_FILES_TO_WRITE
    # global ERROR_FILES_FULL_PATH
    with open(LOG_FILE, 'r') as lf:
        for line in lf:
            if 'ERROR : gene' in line:
                before_keyword, keyword, after_keyword = line.partition('ERROR : gene ')
                gene_name = after_keyword.replace('\n', '')
                ERROR_FILES_TO_WRITE.add(gene_name)
                logging.info("adding error for gene {} from log file {}".format(LOG_FILE, gene_name))
            elif 'OK for file ' in line:
                before_keyword, keyword, after_keyword = line.partition('OK for file ')
                gene_name = after_keyword.replace('\n', '')
                CORRECT_FILES_TO_WRITE.add(gene_name)
                logging.info("adding correct for gene {} from log file {}".format(LOG_FILE, gene_name))


def write_correct_and_error_files(result_file):
    global CORRECT_FILES_TO_WRITE
    global ERROR_FILES_TO_WRITE
    logging.info("CORRECT_FILES_TO_WRITE {}".format(len(CORRECT_FILES_TO_WRITE)))
    logging.info("ERROR_FILES_TO_WRITE {}".format(len(ERROR_FILES_TO_WRITE)))
    writer = pd.ExcelWriter(result_file, engine='openpyxl')
    df_corr = pd.DataFrame({'Gene name': list(CORRECT_FILES_TO_WRITE)})
    df_corr.to_excel(writer, sheet_name='correct files', index=False)
    df_err = pd.DataFrame({'Gene name': list(ERROR_FILES_TO_WRITE)})
    df_err.to_excel(writer, sheet_name='exception files', index=False)
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for Gblocks'
                                    'FASTA formats are accepted', nargs='?', required=True)
    parser.add_argument('--e', help='Input FASTA file extension [fasta, fna, ffn, faa, frn, fa, fas]', nargs='?',
                        default='fas')
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
        pool = multiprocessing.Pool(processes=threads)
        in_dir = args.i
        ext = args.e
        executable_path = args.exec
        input_tuples = list(parse_dir(in_dir, ext))
        len_inputs = len(input_tuples)
        logger.info("Path to the folder with input files for Gblocks: {}\nExecutable path: {}".
                    format(in_dir, executable_path))
        i = pool.starmap_async(launch_gblocks, zip(input_tuples, len_inputs * [args.auto],
                                                   len_inputs * [executable_path], len_inputs * [logger]))
        i.wait()
        i.get()
    except BaseException as e:
        logger.exception("Unexpected error: {}".format(e))
        raise e

    get_correct_and_errors_from_log_file()
    resulting_file = os.path.join(in_dir, 'gblocks_summary.xlsx')
    write_correct_and_error_files(resulting_file)
    # move_error_files(out_dir)
    logging.info("The work has been completed, summary has been written to {}".format(resulting_file))
