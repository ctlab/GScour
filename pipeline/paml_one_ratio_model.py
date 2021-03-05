#!/usr/bin/env python
import argparse
import multiprocessing
import subprocess
from subprocess import TimeoutExpired, SubprocessError
from subprocess import Popen, PIPE
from Bio.Phylo.PAML import codeml
import logging
import os
import re


BROKEN_FILES = list()
PROCESSED_FILES = list()
WRITE_CTL_FILE = 0
LOG_FILE = "paml_one_ratio.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
counter = None


def get_tree_path(trees_folder, species_folder_name):
    for tree in os.scandir(trees_folder):
        if tree.name.split('.')[0] == species_folder_name:
            return tree.name


def get_input_items(folder_in, trees_folder):
    """ parse root folder with files for paml
    parse tree_folder to get appropriate tree """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):
            tree_name = get_tree_path(trees_folder, species_folder.name)
            tree_path = os.path.join(trees_folder, tree_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if infile.split('.')[-1] == 'phy' and infile.split('.')[0].isnumeric():
                            yield folder_in, species_folder.name, item.name, infile, tree_path


def init_indicators(args_counter):
    """store the counter for later use"""
    global counter
    counter = args_counter


def set_one_ratio_model(infile, phylo_tree, personal_dir):
    file_out_path = infile.replace('.phy', '_one_ratio.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylo_tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=1)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=0)
    cml.set_options(NSsites=[0])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=0)
    cml.set_options(omega=.4)
    cml.set_options(fix_alpha=1)
    cml.set_options(alpha=0.)
    cml.set_options(Malpha=0)
    cml.set_options(ncatG=8)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=1)
    cml.set_options(Small_Diff=.5e-6)
    cml.set_options(cleandata=0)
    cml.set_options(method=0)
    cml.write_ctl_file()
    return file_out_path


def run_codeml(input_tuple, exec_path):
    global PROCESSED_FILES
    global counter
    global BROKEN_FILES
    folder_in, species_folder, item_folder, infile, phylogeny_tree_path = input_tuple
    item_folder_path = os.path.join(folder_in, species_folder, item_folder)
    try:
        file_number = (re.search(r"(\d+).phy", infile)).group(1)
    except AttributeError as e:
        logging.info("Please check file .phy {}: cause an error {}".format(item_folder_path, e.args))
        return
    os.chdir(item_folder_path)
    logging.info("Working with {}".format(file_number))
    infile_path = os.path.join(item_folder_path, infile)
    file_out_path = set_one_ratio_model(infile_path, phylogeny_tree_path, item_folder_path)
    if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
        logging.info("Not null size result file {} already exists for file_number {}".format(file_out_path,
                                                                                             file_number))
        return
    p = subprocess.Popen(exec_path, stdin=PIPE, stdout=PIPE)  # '/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml'
    try:
        p.wait(timeout=16000)  # Timeout in seconds: about 4 hours
        if not p.poll():
            raise SubprocessError
        if os.path.getsize(file_out_path) > 0:
            with counter.get_lock():
                counter.value += 1
                logging.info("The work has been done for file {}\nCounter of processed files = {}".
                             format(file_number, counter.value))
            if file_number not in PROCESSED_FILES:
                PROCESSED_FILES.append(file_number)  # TODO: list - to shared variable
                logging.info("PROCESSED_FILES list of length {}: {}".format(len(PROCESSED_FILES), PROCESSED_FILES))
        else:
            logging.info("Null size result file number {}".format(file_number))
            if file_number not in BROKEN_FILES:
                BROKEN_FILES.append(file_number)
                logging.info("BROKEN_FILES list of length {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    except TimeoutExpired as e:
        p.kill()
        logging.info("Killed {}, {}".format(file_number, e))
        if file_number not in BROKEN_FILES:  # TODO: list - to shared variable
            BROKEN_FILES.append(file_number)
            logging.info("BROKEN_FILES of length {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    except SubprocessError:
        logging.info("Not null return code: file number {}".format(file_number))
        if file_number not in BROKEN_FILES:
            BROKEN_FILES.append(file_number)
            logging.info("BROKEN_FILES of length {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))


def main(folder_in, exec_path, trees_folder, threads_number):
    inputs = list(get_input_items(folder_in, trees_folder))  # list of tuples
    len_inputs = len(inputs)
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    counter = multiprocessing.Value('i', 0)
    pool = multiprocessing.Pool(processes=threads_number, initializer=init_indicators,
                                initargs=(counter,))
    i = pool.starmap_async(run_codeml, zip(inputs, len_inputs * [exec_path]))
    i.wait()
    i.get()
    logging.info("Number of files should be analyzed = {}".format(len_inputs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', help='Path to the codeml executable', nargs='?', default="codeml")
    parser.add_argument('--infolder', help='The full path to the folder contains folders with input files for paml',
                        nargs='?')
    parser.add_argument('--tree', help='Path to the folder with trees for paml', nargs='?')
    parser.add_argument('--threads', help='Number of threads to use', nargs='?')
    args = parser.parse_args()
    in_folder = args.infolder
    executable_path = args.e
    tree_folder = args.tree
    threads = int(args.threads)
    logging.info("Path to the folder with input files for paml: {}\nExecutable path: {}\nTree folder: {}\n"
                 "Threads to use = {}".
                 format(in_folder, executable_path, tree_folder, threads))
    try:
        main(in_folder, executable_path, tree_folder, threads)
    except:
        logging.exception("Unexpected error")

    # logging.info("Number of BROCKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    # logging.info("Counter of processed files {}".format(counter.value))
    # logging.info("WRITE_CTL_FILE {}".format(WRITE_CTL_FILE))
    logging.info("The work has been completed")