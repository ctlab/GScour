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
from ctypes import c_wchar_p

BROKEN_FILES = list()
PROCESSED_FILES = list()
WRITE_CTL_FILE = 0
LOG_FILE = "paml_one_ratio.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
counter = None


def parse_dir(folder_in):
    for personal_folder in os.scandir(folder_in):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    yield os.path.join(folder_in, personal_folder, infile)


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


def run_codeml(infile, exec_path, phylogeny_tree):
    global PROCESSED_FILES
    global counter
    global BROKEN_FILES
    personal_dir = os.path.split(infile)[0]
    file_number = (re.search(r"(\d+).phy", infile)).group(1)
    os.chdir(personal_dir)
    logging.info("working with {}".format(file_number))
    file_out_path = set_one_ratio_model(infile, phylogeny_tree, personal_dir)
    if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
        logging.info("Not null size result file {} already exists for file_number {}".format(file_out_path,
                                                                                              file_number))
        return
    p = subprocess.Popen(exec_path, stdin=PIPE, stdout=PIPE)  # '/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml'
    try:
        p.wait(timeout=4000)
        if not p.poll():
            raise SubprocessError
        if os.path.getsize(file_out_path) > 0:
            with counter.get_lock():
                counter.value += 1
                logging.info("The work has been done for file {}\nCounter of processed files = {}".
                             format(file_number, counter.value))
            if file_number not in PROCESSED_FILES:
                PROCESSED_FILES.append(file_number)
                logging.info("PROCESSED_FILES list of length {}: {}".format(len(PROCESSED_FILES), PROCESSED_FILES))
        else:
            logging.info("Null size result file number {}".format(file_number))
            if file_number not in BROKEN_FILES:
                BROKEN_FILES.append(file_number)
                logging.info("BROKEN_FILES list of length {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    except TimeoutExpired as e:
        p.kill()
        logging.info("Killed {}, {}".format(file_number, e))
        if file_number not in BROKEN_FILES:  # to do: list - to shared variable
            BROKEN_FILES.append(file_number)
            logging.info("BROKEN_FILES of length {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    except SubprocessError:
        logging.info("Not null return code, file number {}".format(file_number))
        if file_number not in BROKEN_FILES:
            BROKEN_FILES.append(file_number)
            logging.info("BROKEN_FILES of length {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))


def main(folder_in, phylogeny_tree, exec_path, threads_number):
    """
    for infile in parse_dir(folder_in):
        write_ctl_file(infile, phylogeny_tree)


    for infile in parse_dir(folder_in):
        run_codeml(infile, exec_path)
        """

    inputs = list(parse_dir(folder_in))
    len_inputs = len(inputs)
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    counter = multiprocessing.Value('i', 0)
    # PROCESSED_FILES = multiprocessing.Array('i', range(5900))
    # PROCESSED_FILES = multiprocessing.Array(c_wchar_p, 6000)
    # manager = multiprocessing.Manager()
    # pr = manager.list()
    pool = multiprocessing.Pool(processes=threads_number, initializer=init_indicators,
                                initargs=(counter,))
    i = pool.starmap_async(run_codeml, zip(inputs, len_inputs * [exec_path], len_inputs * [phylogeny_tree]))
    i.wait()
    i.get()
    logging.info("Number of files should be analyzed = {}".format(len_inputs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exec', help='Path to the codeml executable', nargs='?', default="codeml")
    parser.add_argument('--infolder', help='Path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    parser.add_argument('--threads', help='Number of threads to use', nargs='?')
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

    # logging.info("Number of BROCKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    # logging.info("Counter of processed files {}".format(counter.value))
    # logging.info("WRITE_CTL_FILE {}".format(WRITE_CTL_FILE))
    logging.info("The work has been completed")