#!/usr/bin/env python
import argparse
import sys
from Bio.Phylo.PAML import codeml
import logging
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def parse_dir(infolder):
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
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
    cml.set_options(verbose=1)
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
    cml.set_options(verbose=1)
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


def run_paml(infile, tree):
    personal_dir = os.path.split(infile)[0]
    cml, file_out_path = set_null_hypothesis(infile, tree, personal_dir)
    try:
        cml.run(command="codeml", verbose=True) #/home/alina_grf/BIOTOOLS/paml4.9j/bin/
        logging.info("paml out file {} has been written".format(file_out_path))
    except:
        logging.exception("null, sys.exc_info() {}".format(sys.exc_info()))

    cml, file_out_path = set_alternative_hypothesis(infile, tree, personal_dir)
    try:
        cml.run(command="codeml", verbose=True) #/home/alina_grf/BIOTOOLS/paml4.9j/bin/
        logging.info("paml out file {} has been written".format(file_out_path))
    except:
        logging.exception("alter, sys.exc_info() {}".format(sys.exc_info()))


def main(infolder, tree):
    for infile in parse_dir(infolder):
        run_paml(infile, tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='The full path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder, args.tree)
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
