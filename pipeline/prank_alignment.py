#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import logging
import re
import functools, multiprocessing, inspect
import sys
import pandas as pd

CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()
# ERROR_FILES_FULL_PATH = set()
LOG_FILE = "prank_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def parse_initial_dir_create_out_dir(input_dir, extension, output_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == extension:
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


def get_launch_command(infile, final_file_path, outfile_path_without_extension, tree, format_out, aligning):
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile_path_without_extension), 'log'))
    if os.path.isfile(final_file_path):
        raise Exception('final_file_path {} already exists'.format(final_file_path))
    if aligning == 'c':
        aligning = 'codon'
    elif aligning == 't':
        aligning = 'translate'
    if tree:  # -translate (standard code); codon alignment with the option -codon (in -codon case
        # be careful about not multiple of three sequences)
        launch = 'prank -d={0} -o={1} -t={2} -{3} -f={4} > {5}'.format(infile,
                                                                             outfile_path_without_extension,
                                                                             tree, aligning, format_out, log_file)

    else:
        launch = 'prank -d={0} -o={1} -showtree -{2} -f={3} > {4}'.format(infile,
                                                                          outfile_path_without_extension,
                                                                          aligning,
                                                                          format_out,
                                                                          log_file)
    return launch


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except:
            frames = inspect.trace()
            argvalues = inspect.getargvalues(frames[0][0])
            logging.exception("argvalues.locals['args'] {}".format(argvalues.locals['args']))
            file_name = argvalues.locals['args'][0][2]
            gene_name = file_name.split('.')[0]
            logging.exception("gene {}".format(gene_name))
            logging.exception("{} exception in file {}{}/{}".format(sys.exc_info()[0], argvalues.locals['args'][0][0],
                                                                    argvalues.locals['args'][0][1],
                                                                    gene_name))

    return wrapped_func


@trace_unhandled_exceptions
def launch_prank(input_tuple, folder_out, tree, format_out, aligning):
    input_dir, species_folder, infile_name_extended = input_tuple
    infile_path = os.path.join(input_dir, species_folder, infile_name_extended)
    file_gene_name = infile_name_extended.split('.')[0]
    outfile_folder = os.path.join(folder_out, species_folder)
    outfile_path_without_extension = os.path.join(outfile_folder, file_gene_name)
    final_file_name = get_final_file_full_name(outfile_path_without_extension, format_out)

    launch_command = get_launch_command(infile_path, final_file_name,
                                        outfile_path_without_extension, tree,
                                        format_out, aligning)
    os.system(launch_command)


def get_errors_from_log_file():
    global ERROR_FILES_TO_WRITE
    # global ERROR_FILES_FULL_PATH
    with open(LOG_FILE, 'r') as lf:
        for line in lf:
            if 'ERROR : gene' in line:
                before_keyword, keyword, after_keyword = line.partition('ERROR : gene ')
                gene_name = after_keyword.replace('\n', '')
                ERROR_FILES_TO_WRITE.add(gene_name)
                logging.info("adding error for gene {} from log file {}".format(LOG_FILE, gene_name))
            # if 'exception in file ' in line:
            #     before_keyword, keyword, after_keyword = line.partition('exception in file ')
            #     error_file_full_path = after_keyword.replace('\n', '')
            #     ERROR_FILES_FULL_PATH.add(error_file_full_path)


def get_correct_errors_from_prank_log_file(output_dir):
    global CORRECT_FILES_TO_WRITE
    global ERROR_FILES_TO_WRITE
    # global ERROR_FILES_FULL_PATH
    for species_folder in os.scandir(output_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'log':
                    file_gene_name = infile.name.split('.')[0]
                    if os.path.isfile(infile) and os.path.getsize(infile) > 0:
                        with open(infile, 'r') as f_log:
                            result_line = f_log.readlines()[-2]
                            logging.info("result_line {}".format(result_line))
                            if re.search("Analysis done", result_line):
                                logging.info("OK for file".format(infile))
                                CORRECT_FILES_TO_WRITE.add(file_gene_name)
                            else:
                                ERROR_FILES_TO_WRITE.add(file_gene_name)
                                # ERROR_FILES_FULL_PATH.add(os.path.join(output_dir, species_folder.name,
                                # file_gene_name))
                                if file_gene_name in CORRECT_FILES_TO_WRITE:
                                    CORRECT_FILES_TO_WRITE.remove(file_gene_name)
                                logging.error("FAIL for file {}: {}".format(
                                    file_gene_name, result_line))
                    else:
                        ERROR_FILES_TO_WRITE.add(file_gene_name)
                        # ERROR_FILES_FULL_PATH.add(os.path.join(output_dir, species_folder.name, file_gene_name))
                        if file_gene_name in CORRECT_FILES_TO_WRITE:
                            CORRECT_FILES_TO_WRITE.remove(file_gene_name)
                        logging.error("Wrong log file for {}".format(file_gene_name))
                    if file_gene_name in CORRECT_FILES_TO_WRITE and file_gene_name in ERROR_FILES_TO_WRITE:
                        CORRECT_FILES_TO_WRITE.remove(file_gene_name)


def write_correct_and_error_files(output_dir):
    global CORRECT_FILES_TO_WRITE
    global ERROR_FILES_TO_WRITE
    get_errors_from_log_file()
    get_correct_errors_from_prank_log_file(output_dir)
    logging.info("CORRECT_FILES_TO_WRITE {} {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))
    logging.info("ERROR_FILES_TO_WRITE {} {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))
    resulting_file = os.path.join(output_dir, 'prank_summary.xlsx')
    writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
    df_corr = pd.DataFrame({'Gene name': list(CORRECT_FILES_TO_WRITE)})
    df_corr.to_excel(writer, sheet_name='correct files', index=False)
    df_err = pd.DataFrame({'Gene name': list(ERROR_FILES_TO_WRITE)})
    df_err.to_excel(writer, sheet_name='exception files', index=False)
    logging.info("Summary has been written to {}, error files can be moved manually".format(resulting_file))
    writer.save()


# def move_error_files(folder_out): global ERROR_FILES_FULL_PATH logging.info("ERROR_FILES_FULL_PATH {}".format(
# ERROR_FILES_FULL_PATH)) folder_for_errors = os.path.join(folder_out, 'prank_errors') for gene_name_file_path in
# ERROR_FILES_FULL_PATH: nuc_file = gene_name_file_path + ".best.nuc.fas" if os.path.isfile(nuc_file): os.replace(
# os.path.join(nuc_file), os.path.join(folder_for_errors, nuc_file)) pep_file = gene_name_file_path + ".best.pep.fas"
# if os.path.isfile(pep_file): os.replace(os.path.join(folder_out, pep_file), os.path.join(folder_for_errors,
# pep_file)) log_file = gene_name_file_path + ".log" if os.path.isfile(log_file): os.replace(os.path.join(folder_out,
# log_file), os.path.join(folder_for_errors, log_file)) logging.info("nuc {} pep {} log {}".format(nuc_file,
# pep_file, log_file)) # TODO trancfer outfilepath to decorator


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for prank',
                        nargs='?', required=True)
    parser.add_argument('--e', help='Input FASTA file extension [fasta, fna, ffn, faa, frn, fa]', nargs='?',
                        default='fna')
    parser.add_argument('--o', help='Path to the folder with output files of prank, if it does not exist, it will be'
                                    ' created automatically', nargs='?', required=True)
    parser.add_argument('--tree', help='Path to the tree, exclude if there is no tree', nargs='?', default="")
    parser.add_argument('--f', help='Output format', nargs='?',
                        choices=["fasta", "phylipi", "phylips", "paml", "nexus"], default="fasta")
    parser.add_argument('--a', help='Aligning option: \'c\' for codon, \'t\' for translate, ',
                        choices=['t', 'c'], default="t")
    parser.add_argument('--threads', help='Number of threads', nargs='?', required=True)
    args = parser.parse_args()
    threads = int(args.threads)
    in_dir = args.i
    ext = args.e
    out_dir = args.o
    output_format = args.f
    align = args.a
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads)
        input_tuples = list(parse_initial_dir_create_out_dir(in_dir, ext, out_dir))
        len_inputs = len(input_tuples)
        logging.info("files for analysis {}".format(len_inputs))
        logging.info("Options in use:\n Path to the folder with input files for prank: {}\n"
                     "FASTA format for input files: {}\nPath to the folder with output files of prank: {}\n"
                     "tree: {}\noutput format: {}\nAligning option: {}\nthreads: {}"
                     "".format(in_dir, ext, out_dir, args.tree, output_format, align, threads))

        i = pool.starmap_async(launch_prank, zip(input_tuples, len_inputs * [out_dir],
                                                 len_inputs * [args.tree],
                                                 len_inputs * [output_format], len_inputs * [align]))
        i.wait()
        i.get()
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e

    write_correct_and_error_files(out_dir)
    # move_error_files(out_dir)
    logging.info("The work has been completed")
