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


def run_paml(infile, tree):
    file_out_path = infile.replace('phy', 'out')
    personal_dir = os.path.split(file_out_path)[0]

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
    # cml.print_options()
    try:
        cml.run(command="/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml", verbose=True)
        logging.info("paml one ratio model analysis has been done for file {}".format(infile))
        """
        launching of SWAMP masker in separate file due to the need of change python env from 3 to python2
        
        launch_swamp = '/home/alina_grf/BIOTOOLS/SWAMP-master/SWAMP.py' \
                       ' -i {} -b {} -t {} -w {}'.format(personal_dir, branchcodes, threshold, windowsize)
        os.system(launch_swamp)
        logging.info("SWAMP analysis has been done for file {}".format(infile))
        """
    except:
        logging.exception("sys.exc_info() {0}, outfile {1}".format(sys.exc_info(), file_out_path))


def main(infolder, tree):
    for infile in parse_dir(infolder):
        run_paml(infile, tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder, args.tree)
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
