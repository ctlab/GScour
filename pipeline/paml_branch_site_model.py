#!/usr/bin/env python
import argparse
import re
import sys
from Bio.Phylo.PAML import codeml
import os
import logging
from subprocess import STDOUT, check_output


WRITTEN_FILES = 0
EXCEPTION_NUMBER = 0
FILES_NUMBER = 0
BROCKEN_FILES_NULL = list()
BROCKEN_FILES_ALTER = list()
LOG_FILE = "paml_branch_site.log"
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

"""


def parse_dir(infolder):
    global FILES_NUMBER
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    FILES_NUMBER += 1
                    yield os.path.join(infolder, personal_folder, infile)


def set_alternative_hypothesis(infile, tree, personal_dir):
    file_out_path = infile.replace('.phy', '_alter1.out')
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
    return cml, file_out_path


def set_null_hypothesis(infile, tree, personal_dir):
    file_out_path = infile.replace('.phy', '_null1.out')
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
    cml.set_options(fix_omega=1)
    cml.set_options(omega=1)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    return cml, file_out_path


def run_paml(infile, tree, executable_path):
    global WRITTEN_FILES, EXCEPTION_NUMBER
    file_number = (re.search(r"(\d+).phy", infile)).group(1)
    personal_dir = os.path.split(infile)[0]
    cml, file_out_path = set_null_hypothesis(infile, tree, personal_dir)
    try:
        output = check_output(cml.run(command=executable_path, verbose=True),
                              stderr=STDOUT, timeout=120)
        logging.info("paml out file {} has been written".format(file_out_path))
        WRITTEN_FILES += 1
    except:
        logging.exception("null, infile {}, sys.exc_info() {}".format(infile, sys.exc_info()))
        EXCEPTION_NUMBER += 1
        if file_number not in BROCKEN_FILES_NULL:
            BROCKEN_FILES_NULL.append(file_number)

    cml, file_out_path = set_alternative_hypothesis(infile, tree, personal_dir)
    try:
        output = check_output(cml.run(command="/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml", verbose=True),
                              stderr=STDOUT, timeout=120)
        logging.info("paml out file {} has been written".format(file_out_path))
        WRITTEN_FILES += 1
    except:
        logging.exception("alter, infile {},  sys.exc_info() {}".format(infile, sys.exc_info()))
        EXCEPTION_NUMBER += 1
        if file_number not in BROCKEN_FILES_ALTER:
            BROCKEN_FILES_ALTER.append(file_number)


def main(infolder, tree, executable_path):
    for infile in parse_dir(infolder):
        run_paml(infile, tree, executable_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exec', help='Path to the codeml executable', nargs='?', default="codeml")
    parser.add_argument('--infolder', help='The full path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    args = parser.parse_args()
    infolder = args.infolder
    executable_path = args.exec
    tree = args.tree
    logging.info("Path to the folder with input files for paml: {}\nPath to the tree: {}\nExecutable path: {}".
                 format(infolder, tree, executable_path))
    try:
        main(infolder, tree, executable_path)
        logging.info("Number of files have been analyzed: {}".format(FILES_NUMBER))
        logging.info("Number of written files: {}".format(WRITTEN_FILES))
        logging.info("Number of exceptions: {}".format(EXCEPTION_NUMBER))
    except:
        logging.exception("Unexpected error")
        if BROCKEN_FILES_NULL:
            logging.warning("BROCKEN_FILES_NULL: {}".format(BROCKEN_FILES_NULL))
        if BROCKEN_FILES_ALTER:
            logging.warning("BROCKEN_FILES_ALTER: {}".format(BROCKEN_FILES_ALTER))

    logging.info("The work has been completed")
