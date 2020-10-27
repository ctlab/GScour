#!/usr/bin/env python
import argparse
import subprocess
from subprocess import check_output, STDOUT, TimeoutExpired
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


def parse_dir(infolder):
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    yield os.path.join(infolder, personal_folder, infile)


def set_one_ratio_model(infile, tree, personal_dir):
    file_out_path = infile.replace('.phy', '_one_ratio.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=tree,
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
    #cml.write_ctl_file()
    return cml, file_out_path


def write_ctl_file(infile, tree):
    global PROCESSED_FILES
    global WRITE_CTL_FILE
    personal_dir = os.path.split(infile)[0]
    cml, file_out_path = set_one_ratio_model(infile, tree, personal_dir)
    os.chdir(personal_dir)
    cml.write_ctl_file()
    WRITE_CTL_FILE += 1


def run_codeml(infile):
    personal_dir = os.path.split(infile)[0]
    file_number = (re.search(r"(\d+).phy", infile)).group(1)
    os.chdir(personal_dir)
    logging.info("working with {}".format(file_number))
    p = subprocess.Popen('/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml', stdin=PIPE, stdout=PIPE)
    try:
        p.wait(timeout=20)
        logging.info("The work has been done for file {}".format(file_number))
        if file_number not in PROCESSED_FILES:
            PROCESSED_FILES.append(file_number)
    except TimeoutExpired as e:
        p.kill()
        logging.info("Killed {}".format(file_number))
        if file_number not in BROKEN_FILES:
            BROKEN_FILES.append(file_number)


def main(infolder, tree):
    for infile in parse_dir(infolder):
        write_ctl_file(infile, tree)
    for infile in parse_dir(infolder):
        run_codeml(infile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder, args.tree)
    except:
        logging.exception("Unexpected error")
        logging.warning("NUMBER OF BROCKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
        logging.info("NUMBER OF PROCESSED FILES {}:{}".format(len(PROCESSED_FILES), PROCESSED_FILES))

    logging.warning("NUMBER OF BROCKEN_FILES {}: {}".format(len(BROKEN_FILES), BROKEN_FILES))
    logging.info("NUMBER OF PROCESSED FILES {}:{}".format(len(PROCESSED_FILES), PROCESSED_FILES))
    logging.info("WRITE_CTL_FILE {}".format(WRITE_CTL_FILE))
    logging.info("The work has been completed")
