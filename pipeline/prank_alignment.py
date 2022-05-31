#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import pandas as pd

CORRECT_FILES = list()
ERROR_FILES = list()
CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()
exception_number = 0
counter = None
LOG_FILE = "prank_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(counter_args):
    """store the counter for later use"""
    global counter
    counter = counter_args


def parse_dir(input_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fna':
                    yield input_dir, species_folder.name, infile.name


def get_final_file_path(outfile_path_without_extension, format_out):
    if format_out in ["phylipi", "phylips", "paml"]:
        final_file_path = '{}.best.phy'.format(outfile_path_without_extension)
    elif format_out == "fasta":
        final_file_path = '{}.best.fas'.format(outfile_path_without_extension)
    else:
        final_file_path = '{}.best.nex'.format(outfile_path_without_extension)
    return final_file_path


def get_launch_command(infile, final_file_path, outfile_path_without_extension, tree, format_out):
    log_file = os.path.join("{}.{}".format(os.path.abspath(outfile_path_without_extension), 'log'))
    if os.path.isfile(final_file_path):
        raise Exception('final_file_path {} already exists'.format(final_file_path))
    if tree and format_out:  # -translate (standard code); codon alignment with the option -codon (in -codon case
        # be careful about not multiple of three sequences)
        launch = 'prank -d={0} -o={1} -t={2} -translate -f={3} > {4}'.format(infile, outfile_path_without_extension,
                                                                             tree, format_out, log_file)
    elif tree:
        launch = 'prank -d={0} -o={1} -t={2} -translate > {3}'.format(infile, outfile_path_without_extension,
                                                                      tree, log_file)
    elif format_out and not tree:
        launch = 'prank -d={0} -o={1} -showtree -translate -f={2} > {3}'.format(infile,
                                                                                outfile_path_without_extension,
                                                                                format_out,
                                                                                log_file)
    else:
        launch = 'prank -d={0} -o={1} -showtree -translate > {2}'.format(infile,
                                                                         outfile_path_without_extension,
                                                                         log_file)
    return launch


def launch_prank(input_tuple, folder_out, tree, format_out):
    global counter
    global exception_number
    input_dir, species_folder, infile_name_extended = input_tuple
    infile_path = os.path.join(input_dir, species_folder, infile_name_extended)
    file_gene_name = infile_name_extended.split('.')[0]
    outfile_path_without_extension = os.path.join(folder_out, file_gene_name)
    final_file_path = get_final_file_path(outfile_path_without_extension, format_out)
    try:
        global CORRECT_FILES
        global ERROR_FILES
        # file_gene_name = re.search(r'(\d+)\.', infile_name_extended).group(1)
        launch_command = get_launch_command(infile_path, final_file_path, outfile_path_without_extension, tree,
                                            format_out)
        if not os.system(launch_command):
            # logging.info("prank completed task for file {}".format(file_gene_name))
            if file_gene_name not in CORRECT_FILES:
                CORRECT_FILES.append(file_gene_name)
                with counter.get_lock():
                    counter.value += 1  # TODO: wrong number

                    # logging.info("Counter (ALIGNED_FILES) = {}".format(counter.value)) # TODO: multiprocessing
    except BaseException as err:
        ERROR_FILES.append(file_gene_name)
        logging.exception("Infile {}, - Unexpected error: {}".format(infile_name_extended, err))
        exception_number += 1
        logging.info("exception number = {}".format(exception_number))
    return CORRECT_FILES, ERROR_FILES


def gather_correct(correct):
    logging.info("res {} {}".format(type(correct), correct))
    for tup in correct:
        for corr in tup:
            for i in corr:
                CORRECT_FILES_TO_WRITE.add(i)
    logging.info("CORRECT_FILES_TO_WRITE {} {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))


def gather_errors(exception):
    logging.info("exception {} {}".format(type(exception), exception))
    for tup in exception:
        for e in tup:
            for i in e:
                ERROR_FILES_TO_WRITE.add(i)
    logging.info("ERROR_FILES_TO_WRITE {} {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))


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
    global exception_list
    exception_list = list()
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        in_dir = args.i
        input_tuples = list(parse_dir(in_dir))
        len_inputs = len(input_tuples)
        out_dir = args.o
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        output_format = args.f
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
        # corrects, errors = returned_tuple
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    # resulting_file = os.path.join(out_dir, 'prank_summary.xlsx')
    # logging.info("res file {}".format(resulting_file))
    # writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
    # df_corr = pd.DataFrame({'Gene name': CORRECT_FILES})
    # df_corr.to_excel(writer, sheet_name='correct files', index=False)
    # df_err = pd.DataFrame({'Gene name': ERROR_FILES})
    # df_corr.to_excel(writer, sheet_name='error files', index=False)
    # writer.save()
    logging.info("The work has been completed")
