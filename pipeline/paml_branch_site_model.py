#!/usr/bin/env python
import argparse
import multiprocessing
import re
import subprocess
import traceback
from subprocess import SubprocessError, TimeoutExpired
from Bio.Phylo.PAML import codeml
import os
import logging

BROCKEN_FILES_NULL = list()
BROCKEN_FILES_ALTER = list()
excep_null_counter = None
excep_alter_counter = None
PROCESSED_FILES_NULL = list()
PROCESSED_FILES_ALTER = list()
processed_null_counter = None
processed_alter_counter = None

LOG_FILE = "paml_branch_site.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)

"""There are two hypothesis:
H0 (The null model for Branch-site model A): 
    Model A1: model = 2, NSsites = 2, fix_omega = 1, omega = 1
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
H1 (Alternative model, Model A: model = 2, NSsites = 2, fix_omega = 0 ): 
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
    
From readme of paml example lysozymeLarge.ctl:
Alternative hypothesis (branch site model A, with w2 estimated):
model = 2    NSsites = 2   fix_omega = 0   omega = 1.5 (or any value > 1)

As the branch-site model is
known to cause computational difficulties for the numer-
ical optimization algorithm in PAML, each analysis is con-
ducted three times with different starting values to ensure
that the global peak is found (Statistical Properties of the Branch-Site Test of Positive
Selection, Ziheng Yang and Mario dos Reis)
"""


def get_tree_path(trees_folder, species_folder_name):
    for tree in os.scandir(trees_folder):
        if tree.name.split('.')[0] == species_folder_name:
            return tree.name


def get_input_items(folder_in, trees_folder):
    """ parse root folder with files for paml
    parse tree_folder to get appropriate tree """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):
            logging.info("working with species folder {}".format(species_folder.name))
            tree_name = get_tree_path(trees_folder, species_folder.name)
            tree_path = os.path.join(trees_folder, tree_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if infile.split('.')[-1] == 'phy' and not re.search('[a-zA-Z]', infile.split('.')[0]):
                            yield folder_in, species_folder.name, item.name, infile, tree_path


def init_indicators(null_args, alter_args, excep_null, excep_alt):
    """store the counter for later use"""
    global processed_null_counter, processed_alter_counter
    global excep_null_counter, excep_alter_counter
    processed_null_counter = null_args
    processed_alter_counter = alter_args
    excep_null_counter = excep_null
    excep_alter_counter = excep_alt


def set_alternative_hypothesis(infile, phylo_tree, personal_dir):
    file_out_path = infile.replace('.phy', '_alter1.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylo_tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=0)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=2)
    cml.set_options(NSsites=[2])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=0)
    cml.set_options(omega=2)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    cml.write_ctl_file()
    return cml, file_out_path


def set_null_hypothesis(infile, phylo_tree, personal_dir):
    file_out_path = infile.replace('.phy', '_null1.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylo_tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=0)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=2)
    cml.set_options(NSsites=[2])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=1)
    cml.set_options(omega=1)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    cml.write_ctl_file()
    return cml, file_out_path


def run_paml(input_tuple, exec_path, hypothesis_type, overwrite_flag):
    folder_in, species_folder, item_folder, infile, phylogeny_tree_path = input_tuple
    item_folder_path = os.path.join(folder_in, species_folder, item_folder)
    infile_path = os.path.join(item_folder_path, infile)
    try:
        file_number = (re.search(r"(\d+).phy", infile)).group(1)
    except AttributeError as err:
        logging.info("There is no .phy for {}:\n{}".format(infile_path, err.args))
        return
    os.chdir(item_folder_path)
    logging.info("Working with {}".format(file_number))

    if hypothesis_type == "null":
        global BROCKEN_FILES_NULL, PROCESSED_FILES_NULL
        global excep_null_counter, processed_null_counter
        broken_files = BROCKEN_FILES_NULL
        processed_files = PROCESSED_FILES_NULL
        excep_counter = excep_null_counter
        processed_counter = processed_null_counter
        cml, file_out_path = set_null_hypothesis(infile_path, phylogeny_tree_path, item_folder_path)
    elif hypothesis_type == "alter":
        global BROCKEN_FILES_ALTER, PROCESSED_FILES_ALTER
        global excep_alter_counter, processed_alter_counter
        broken_files = BROCKEN_FILES_ALTER
        processed_files = PROCESSED_FILES_ALTER
        excep_counter = excep_alter_counter
        processed_counter = processed_alter_counter
        cml, file_out_path = set_alternative_hypothesis(infile_path, phylogeny_tree_path, item_folder_path)
    else:
        logging.warning("Check the type of hypothesis: null, alter")
        return

    if not overwrite_flag:
        if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
            with open(file_out_path, 'r') as o_f:
                if "Time used" in o_f.readlines()[-1]:
                    logging.info("The work has already been done for file {}...continue...".format(file_out_path))
                    return

    p = subprocess.Popen(exec_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        return_code = p.wait(timeout=500)  # Timeout in seconds
        stdout, stderr = p.communicate()
        logging.info("Return code={} for file number {}, stderr: {}".format(return_code, file_number, stderr))
        if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
            with open(file_out_path, 'r') as o_f:
                if "Time used" in o_f.readlines()[-1]:
                    with processed_counter.get_lock():
                        processed_counter.value += 1
                    logging.info("The work has been done for file {}\nCounter of processed files = {}".
                                 format(file_number, processed_counter.value))
                    if file_number not in processed_files:
                        processed_files.append(file_number)  # TODO: list - to shared variable
                        # logging.info(
                        #     "PROCESSED_FILES list of length {}: {}".format(len(processed_files), processed_files))
                else:
                    logging.warning("The work has not been finished for file number {}".format(file_number))
                    if file_number not in broken_files:
                        broken_files.append(file_number)
                        # logging.info("BROKEN_FILES list of length {}: {}".format(len(broken_files), broken_files))
    except TimeoutExpired as err:
        p.kill()
        with excep_counter.get_lock():
            excep_counter.value += 1
        file_id = "{}/{}".format(species_folder, item_folder)
        logging.warning("Killed {} hypothesis_type {}, {}, {},\nException_counter={}".format(file_id, hypothesis_type,
                                                                                          err.args, err,
                                                                                          excep_counter.value))
        if file_id not in broken_files:  # TODO: list - to shared variable
            broken_files.append(file_id)
            # logging.info("BROKEN_FILES of length {}: {}".format(len(broken_files), broken_files))


def main(folder_in, trees_folder, exec_path, number_of_threads, overwrite_flag):
    inputs = list(get_input_items(folder_in, trees_folder))
    len_inputs = len(inputs)
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)

    null_processed_counter = multiprocessing.Value('i', 0)
    alter_processed_counter = multiprocessing.Value('i', 0)
    null_exception_counter = multiprocessing.Value('i', 0)
    alter_exception_counter = multiprocessing.Value('i', 0)
    # PROCESSED_FILES = multiprocessing.Array('i', range(5900))
    # PROCESSED_FILES = multiprocessing.Array(c_wchar_p, 6000)
    # manager = multiprocessing.Manager()
    # pr = manager.list()
    pool = multiprocessing.Pool(processes=number_of_threads, initializer=init_indicators,
                                initargs=(null_processed_counter, alter_processed_counter, null_exception_counter,
                                          alter_exception_counter,))
    i = pool.starmap_async(run_paml, zip(inputs, len_inputs * [exec_path],
                                         len_inputs * ["null"], len_inputs * [overwrite_flag]))
    i.wait()
    i.get()
    i = pool.starmap_async(run_paml, zip(inputs, len_inputs * [exec_path],
                                         len_inputs * ["alter"], len_inputs * [overwrite_flag]))
    i.wait()
    i.get()
    logging.info("Number of files should be analyzed: {}".format(len_inputs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', help='Path to the codeml executable', nargs='?', default="codeml")
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?')
    parser.add_argument('--tree', help='Path to the folder with trees for paml', nargs='?')
    parser.add_argument('--threads', help='Number of threads to use', nargs='?')
    parser.add_argument('--rework', help='"y" if overwrite existing files, default "n"', nargs='?', default='n')
    args = parser.parse_args()
    in_folder = args.i
    executable_path = args.e
    tree_folder = args.tree
    threads = int(args.threads)
    if args.rework == 'y':
        rework = True
    else:
        rework = False
    logging.info("Path to the folder with input files for paml: {}\nPath to the tree: {}\nExecutable path: {}\n"
                 "Threads to use = {}, rework = {}".
                 format(in_folder, tree_folder, executable_path, threads, rework))
    try:
        main(in_folder, tree_folder, executable_path, threads, rework)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        # logging.info("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))
        # if BROCKEN_FILES_NULL:
        #     logging.warning("BROCKEN_FILES_NULL: {}".format(BROCKEN_FILES_NULL))
        # if BROCKEN_FILES_ALTER:
        #     logging.warning("BROCKEN_FILES_ALTER: {}".format(BROCKEN_FILES_ALTER))
    # logging.info("Counter of processed files {}".format(counter.value))
    # logging.info("Number of exceptions: {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")
