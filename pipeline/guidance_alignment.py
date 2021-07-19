#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import multiprocessing
import os
import logging
import re
import shutil

PROCESSED_FILES = list()
EXCEPTION_NUMBER = 0
counter = None
LOG_FILE = "guidance_alignment.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def init_counter(args_counter):
    """store the counter for later use"""
    global counter
    counter = args_counter


def parse_dir(input_dir):
    for species_folder in os.scandir(input_dir):
        if os.path.isdir(species_folder):
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fna':
                    yield input_dir, species_folder.name, infile.name


def find_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def launch_guidance(input_tuple, folder_out, number_of_threads, executable_path):
    global counter
    global PROCESSED_FILES
    global EXCEPTION_NUMBER
    input_dir, species_folder_name, infile_name = input_tuple
    infile_path = os.path.join(input_dir, species_folder_name, infile_name)
    file_number = re.search(r'(\d+)\.', infile_name).group(1)
    personal_dir_full_out = os.path.join(folder_out, species_folder_name, file_number)
    if not os.path.isdir(personal_dir_full_out):
        os.makedirs(personal_dir_full_out)
    if not os.path.isdir(os.path.join(folder_out, "cleansed", species_folder_name)):
        os.makedirs(os.path.join(folder_out, "cleansed", species_folder_name))
    cleansed_file_path = os.path.join(folder_out, "cleansed", species_folder_name, "{}.{}".format(file_number, 'fna'))
    """
    Required parameters:
    --seqFile: Input sequence file in FASTA format
    --msaProgram: Which MSA program to use
    --seqType: Type of sequences for alignment (amino acids, nucleotides, or codons)
    --outDir: Output directory that will be created automatically and hold all output files
     [please provide full (and not relative) path]
    --bootstraps: Number of bootstrap iterations (only for GUIDANCE). Default=100
    --prank: path to prank executable. Default=prank
    """
    msaProgram = 'PRANK'
    launch = 'perl {} --seqFile {} --msaProgram {} ' \
             '--seqType nuc ' \
             '--proc_num {} --outDir {} --bootstraps 30 --outOrder as_input'.format(executable_path, infile_path,
                                                                                    msaProgram,
                                                                                    number_of_threads,
                                                                                    personal_dir_full_out)
    print("launch", launch)
    """
    final result will be recorded in file MSA.{msaProgram_name}.Without_low_SP_Col.With_Names
    """

    os.system(launch)
    logging.info("Guidance completed task for file {}".format(file_number))
    if file_number not in PROCESSED_FILES:
        PROCESSED_FILES.append(file_number)  # TODO: to shared variables
        with counter.get_lock():
            counter.value += 1
            logging.info("Counter (PROCESSED_FILES) = {}\nList of PROCESSED_FILES: {}".
                         format(counter.value, PROCESSED_FILES))
    resulting_file_name = 'MSA.{}.Without_low_SP_Col.With_Names'.format(msaProgram)
    resulting_file_path = find_file(resulting_file_name, personal_dir_full_out)
    if resulting_file_path:
        shutil.copy(resulting_file_path, cleansed_file_path)
        logging.info("File {} recorded".format(cleansed_file_path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files (.fna) for guidance',
                        nargs='?')
    parser.add_argument('--o', help='Path to the folder with output files of guidance', nargs='?')
    parser.add_argument('--exec', help='Path to the guidance executable "guidance.pl"', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    threads = int(args.threads)
    in_dir = args.i
    out_dir = args.o
    logging.info("Path to the folder with input files for guidance: {}\n"
                 "Path to the folder with output files of guidance: {}\n".format(in_dir, out_dir))
    try:
        counter = multiprocessing.Value('i', 0)
        multiprocessing.log_to_stderr()
        logger = multiprocessing.get_logger()
        logger.setLevel(logging.INFO)
        pool = multiprocessing.Pool(processes=threads, initializer=init_counter, initargs=(counter,))
        input_tuples = list(parse_dir(in_dir))
        len_inputs = len(input_tuples)
        i = pool.starmap_async(launch_guidance, zip(input_tuples, len_inputs * [out_dir], len_inputs * [threads],
                                                    len_inputs * [args.exec]))
        i.wait()
        i.get()
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))

    logging.info("Number of ALIGNED_FILES = {}".format(counter.value))
    logging.info("Number of exceptions = {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")