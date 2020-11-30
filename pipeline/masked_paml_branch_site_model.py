#!/usr/bin/env python
import argparse
import multiprocessing
import re
import subprocess
from subprocess import TimeoutExpired, SubprocessError
import sys
from Bio.Phylo.PAML import codeml
import logging
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

LOG_FILE = "paml_branch_site_masked.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)

"""There are two hypothesis:
H0: model = 2, NSsites = 2 (branch-site model),
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
H1: model = 2, NSsites = 2 (branch-site model),
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
    
in analyse_pamlout.py:
For analysis performs: ln0, np0 from H0; ln1, np1 from H1; 
                       ΔLRT = 2×(lnL1 - lnL0)
                       n = np1 - np0
                       p_val = 1 - stats.chi2.cdf(ΔLRT, n)
                       
Launch this script for files masked with SWAMP:
1. paml, one-ratio model
2. SWAMP
3. paml, branch-site model

"""


def parse_dir(folder_in):
    for personal_folder in os.scandir(folder_in):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.endswith("_masked.phy"):
                    yield os.path.join(folder_in, personal_folder, infile)


def init_indicators(null_args, alter_args, excep_null, excep_alt):
    """store the counter for later use"""
    global processed_null_counter, processed_alter_counter
    global excep_null_counter, excep_alter_counter
    counter_null = null_args
    counter_alter = alter_args
    excep_null_counter = excep_null
    excep_alter_counter = excep_alt


def set_alternative_hypothesis(infile, tree, personal_dir):
    file_out_path = infile.replace('_masked.phy', '_alter1_masked.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=tree,
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
    cml.set_options(omega=1)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    cml.write_ctl_file()
    return cml, file_out_path


def set_null_hypothesis(infile, phylogeny_tree, personal_dir):
    file_out_path = infile.replace('_masked.phy', '_null1_masked.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylogeny_tree,
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


def run_paml(infile, phylo_tree, exec_path, hypothesis_type):
    personal_dir = os.path.split(infile)[0]
    file_number = (re.search(r"(\d+).phy", infile)).group(1)
    os.chdir(personal_dir)
    logging.info("working with {}".format(file_number))  # '/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml'

    if hypothesis_type == "null":
        global BROCKEN_FILES_NULL, PROCESSED_FILES_NULL
        global excep_null_counter, processed_null_counter
        BROKEN_FILES = BROCKEN_FILES_NULL
        PROCESSED_FILES = PROCESSED_FILES_NULL
        excep_counter = excep_null_counter
        processed_counter = processed_null_counter
        cml, file_out_path = set_null_hypothesis(infile, phylo_tree, personal_dir)
    elif hypothesis_type == "alter":
        global BROCKEN_FILES_ALTER, PROCESSED_FILES_ALTER
        global excep_alter_counter, processed_alter_counter
        BROKEN_FILES = BROCKEN_FILES_ALTER
        PROCESSED_FILES = PROCESSED_FILES_ALTER
        excep_counter = excep_alter_counter
        processed_counter = processed_alter_counter
        cml, file_out_path = set_alternative_hypothesis(infile, phylo_tree, personal_dir)

    if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
        logging.info("{}: Not null size result file {} already exists for file_number {}".
                     format(hypothesis_type, file_out_path, file_number))
        return
    p = subprocess.Popen(exec_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        p.wait(timeout=4000)
        if not p.poll():
            raise SubprocessError
        if os.path.getsize(file_out_path) > 0:
            with processed_counter.get_lock():
                processed_counter.value += 1
                logging.info("{}: The work has been done for file {}\nCounter of processed files = {}".
                             format(hypothesis_type, file_number, processed_counter.value))
            if file_number not in PROCESSED_FILES:
                PROCESSED_FILES.append(file_number)
                logging.info("{}: PROCESSED_FILES list of length {}: {}".format(hypothesis_type, len(PROCESSED_FILES),
                                                                                PROCESSED_FILES))
        else:
            logging.info("{}: Null size result file number {}".format(hypothesis_type, file_number))
            if file_number not in BROKEN_FILES:
                BROKEN_FILES.append(file_number)
                logging.info("{}: BROKEN_FILES list of length {}: {}".format(hypothesis_type, len(BROKEN_FILES),
                                                                             BROKEN_FILES))
    except TimeoutExpired as e:
        p.kill()
        logging.info("{}: Killed {}, {}".format(hypothesis_type, file_number, e))
        with excep_counter.get_lock():
            excep_counter.value += 1
            logging.info("{}: exception counter {}".
                         format(hypothesis_type, excep_counter.value))
        if file_number not in BROKEN_FILES:  # to do: list - to shared variable
            BROKEN_FILES.append(file_number)
            logging.info("{}: BROKEN_FILES list of length {}: {}".format(hypothesis_type, len(BROKEN_FILES),
                                                                         BROKEN_FILES))
    except SubprocessError:
        logging.info("{}: Not null return code, file number {}".format(hypothesis_type, file_number))
        with excep_counter.get_lock():
            excep_counter.value += 1
            logging.info("{}: exception counter {}".
                         format(hypothesis_type, excep_counter.value))
        if file_number not in BROKEN_FILES:
            BROKEN_FILES.append(file_number)
            logging.info("{}: BROKEN_FILES of length {}: {}".format(hypothesis_type, len(BROKEN_FILES), BROKEN_FILES))
"""
    cml, file_out_path = set_alternative_hypothesis(infile, phylo_tree, personal_dir)
    if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
        logging.info("Not null size result file {} already exists for file_number {}".format(file_out_path,
                                                                                             file_number))
        return
    p = subprocess.Popen(exec_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        p.wait(timeout=4000)
        if not p.poll():
            raise SubprocessError
        logging.info("paml out file {} has been written".format(file_out_path))
        with counter.get_lock():
            counter.value += 1
            logging.info("Counter of processed files = {}".format(counter.value))
        if file_number not in PROCESSED_FILES:
            PROCESSED_FILES.append(file_number)
    except subprocess.TimeoutExpired as e:
        p.kill()
        EXCEPTION_NUMBER += 1
        logging.exception("Alternative hypothesis, Killed {}, {}".format(file_number, e))
        if file_number not in BROCKEN_FILES_ALTER:  # to do: list - to shared variable
            BROCKEN_FILES_ALTER.append(file_number)
    except SubprocessError:
        logging.info("Not null return code, file number {}".format(file_number))
        if file_number not in BROCKEN_FILES_NULL:
            BROCKEN_FILES_ALTER.append(file_number)
"""


def main(folder_in, phylogeny_tree, exec_path, number_of_threads):
    """
    for infile in parse_dir(folder_in):
        run_paml(infile, phylogeny_tree, exec_path) """
    inputs = list(parse_dir(folder_in))
    len_inputs = len(inputs)
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    null_processed_counter = multiprocessing.Value('i', 0)
    alter_processed_counter = multiprocessing.Value('i', 0)
    null_excep_counter = multiprocessing.Value('i', 0)
    alter_excep_counter = multiprocessing.Value('i', 0)
    # PROCESSED_FILES = multiprocessing.Array('i', range(5900))
    # PROCESSED_FILES = multiprocessing.Array(c_wchar_p, 6000)
    # manager = multiprocessing.Manager()
    # pr = manager.list()
    pool = multiprocessing.Pool(processes=number_of_threads, initializer=init_indicators,
                                initargs=(null_processed_counter,alter_processed_counter,null_excep_counter,
                                          alter_excep_counter,))
    i = pool.starmap_async(run_paml, zip(inputs, len_inputs * [phylogeny_tree], len_inputs * [exec_path],
                                         len_inputs * ["null"]))
    i.wait()
    i.get()
    i = pool.starmap_async(run_paml, zip(inputs, len_inputs * [phylogeny_tree], len_inputs * [exec_path],
                                         len_inputs * ["alter"]))
    i.wait()
    i.get()
    logging.info("Number of files should be analyzed: {}".format(len_inputs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exec', help='Path to the codeml executable', nargs='?', default="codeml")
    parser.add_argument('--infolder', help='The full path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    parser.add_argument('--threads', help='Number of threads', nargs='?')
    args = parser.parse_args()
    infolder = args.infolder
    executable_path = args.exec
    tree = args.tree
    threads = int(args.threads)
    logging.info("Path to the folder with input files for paml: {}\nPath to the tree: {}\nExecutable path: {}\n"
                 "Threads to use = {}".
                 format(infolder, tree, executable_path, threads))
    try:
        main(infolder, tree, executable_path, threads)
    except:
        logging.exception("Unexpected error")
        if BROCKEN_FILES_NULL:
            logging.warning("BROCKEN_FILES_NULL: {}".format(BROCKEN_FILES_NULL))
        if BROCKEN_FILES_ALTER:
            logging.warning("BROCKEN_FILES_ALTER: {}".format(BROCKEN_FILES_ALTER))
    # logging.info("Counter of processed files {}".format(counter.value))
    # logging.info("Number of exceptions: {}".format(EXCEPTION_NUMBER))
    logging.info("The work has been completed")
