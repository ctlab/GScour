#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import pandas as pd


CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()
CORRECT_FILES = list()
ERROR_FILES = list()
LOG_FILE = "prank_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(counter_args):
    """store the counter for later use"""
    global counter
    counter = counter_args


def parse_initial_dir_create_out_dir(input_dir, output_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fna':
                    if not os.path.isdir(os.path.join(output_dir, species_folder.name)):
                        os.makedirs(os.path.join(output_dir, species_folder.name))
                    yield input_dir, species_folder.name, infile.name


def get_final_file_full_name(outfile_path_without_extension, format_out):
    if format_out in ["phylipi", "phylips", "paml"]:
        final_file_full_name = '{}.best.phy'.format(outfile_path_without_extension)
    elif format_out == "fasta":
        final_file_full_name = '{}.best.fas'.format(outfile_path_without_extension)
    else:
        final_file_full_name = '{}.best.nex'.format(outfile_path_without_extension)
    return final_file_full_name


def get_launch_command(infile, final_file_path, outfile_path_without_extension, tree, format_out):
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile_path_without_extension), 'log'))
    if os.path.isfile(final_file_path):
        raise Exception('final_file_path {} already exists'.format(final_file_path))
    if tree:  # -translate (standard code); codon alignment with the option -codon (in -codon case
        # be careful about not multiple of three sequences)
        launch = 'prank -d={0} -o={1} -t={2} -once -translate -f={3} > {4}'.format(infile,
                                                                                   outfile_path_without_extension,
                                                                                   tree, format_out, log_file)

    else:
        launch = 'prank -d={0} -o={1} -showtree -translate -f={2} > {3}'.format(infile,
                                                                                outfile_path_without_extension,
                                                                                format_out,
                                                                                log_file)
    return launch


def launch_prank(input_tuple, folder_out, tree, format_out):
    input_dir, species_folder, infile_name_extended = input_tuple
    infile_path = os.path.join(input_dir, species_folder, infile_name_extended)
    file_gene_name = infile_name_extended.split('.')[0]
    outfile_folder = os.path.join(folder_out, species_folder)
    outfile_path_without_extension = os.path.join(outfile_folder, file_gene_name)
    final_file_name = get_final_file_full_name(outfile_path_without_extension, format_out)
    try:
        global CORRECT_FILES
        global ERROR_FILES
        launch_command = get_launch_command(infile_path, final_file_name,
                                            outfile_path_without_extension, tree,
                                            format_out)
        if not os.system(launch_command):
            if file_gene_name not in CORRECT_FILES:
                CORRECT_FILES.append(file_gene_name)

    except BaseException as err:
        ERROR_FILES.append(file_gene_name)
        logging.exception("Infile {}, - Unexpected error: {}".format(infile_name_extended, err))
    return CORRECT_FILES, ERROR_FILES


def gather_correct(correct):
    global CORRECT_FILES_TO_WRITE
    for tup in correct:
        if not tup[0]:
            logging.warning("No correct files, check prank log file")
        for correct_item in tup[0]:
            CORRECT_FILES_TO_WRITE.add(correct_item)
    logging.info("CORRECT_FILES_TO_WRITE {} {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))


def gather_errors(exception):
    global ERROR_FILES_TO_WRITE
    for tup in exception:
        for error_item in tup[1]:
            ERROR_FILES_TO_WRITE.add(error_item)
    logging.info("ERROR_FILES_TO_WRITE {} {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))


def write_correct_error_files(output_dir):
    global CORRECT_FILES_TO_WRITE
    global ERROR_FILES_TO_WRITE
    for species_folder in os.scandir(output_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'log':
                    if os.path.isfile(infile) and os.path.getsize(infile) > 0:
                        file_gene_name = infile.name.split('.')[0]
                        with open(infile, 'r') as f_log:
                            result_line = f_log.readlines()[-2]
                            logging.info("result_line {}".format(result_line))
                            if re.search("Analysis done", result_line):
                                logging.info("Log parsing: analysis done for file".format(infile))
                                CORRECT_FILES_TO_WRITE.add(file_gene_name)
                            else:
                                ERROR_FILES_TO_WRITE.add(file_gene_name)
                                if file_gene_name in CORRECT_FILES_TO_WRITE:
                                    CORRECT_FILES_TO_WRITE.remove(file_gene_name)
                                logging.error("Log parsing: the work hasn't been complete for file {}: {}".format(
                                    file_gene_name, result_line))
                    else:
                        ERROR_FILES_TO_WRITE.add(file_gene_name)
                        if file_gene_name in CORRECT_FILES_TO_WRITE:
                            CORRECT_FILES_TO_WRITE.remove(file_gene_name)
                        logging.error("Wrong log file for {}".format(file_gene_name))
    logging.info("CORRECT_FILES_TO_WRITE {} {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))
    logging.info("ERROR_FILES_TO_WRITE {} {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))
    resulting_file = os.path.join(output_dir, 'prank_summary.xlsx')
    logging.info("res file {}".format(resulting_file))
    writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
    df_corr = pd.DataFrame({'Gene name': list(CORRECT_FILES_TO_WRITE)})
    df_corr.to_excel(writer, sheet_name='correct files', index=False)
    df_err = pd.DataFrame({'Gene name': list(ERROR_FILES_TO_WRITE)})
    df_err.to_excel(writer, sheet_name='exception files', index=False)
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files (.fna) for prank',
                        nargs='?', required=True)
    parser.add_argument('--o', help='Path to the folder with output files of prank, if it does not exist, it will be'
                                    ' created automatically', nargs='?', required=True)
    parser.add_argument('--tree', help='Path to the tree, exclude if there is no tree', nargs='?', default="")
    parser.add_argument('--f', help='Output format: ["fasta" (default option for the pipeline, exclude option --f if '
                                    'left by default),'
                                    '"phylipi", "phylips", "paml", "nexus"]', nargs='?', default="fasta")
    parser.add_argument('--threads', help='Number of threads', nargs='?', required=True)
    args = parser.parse_args()
    threads = int(args.threads)
    in_dir = args.i
    out_dir = args.o
    output_format = args.f
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        input_tuples = list(parse_initial_dir_create_out_dir(in_dir, out_dir))
        len_inputs = len(input_tuples)
        if output_format not in ["fasta", "phylipi", "phylips", "paml", "nexus"]:
            raise SyntaxError("Not valid output format, check option --f, -h for help")
        logging.info("Path to the folder with input files for prank: {}\n"
                     "Path to the folder with output files of prank: {}\n"
                     "tree: {}\noutput format: {}\nthreads: {}"
                     "".format(in_dir, out_dir, args.tree, output_format, threads))

        i = pool.starmap_async(launch_prank, zip(input_tuples, len_inputs * [out_dir], len_inputs * [args.tree],
                                                 len_inputs * [output_format]), callback=gather_correct,
                               error_callback=gather_errors)
        i.wait()
        i.get()
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e

    write_correct_error_files(out_dir)
    logging.info("The work has been completed")
