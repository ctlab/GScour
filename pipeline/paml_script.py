import sys

from Bio.Phylo.PAML import codeml
import logging
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
directory_in = '/home/alina_grf/progprojects/data/acfhp_phylip_strict_seqio/'
tree_file = '/home/alina_grf/progprojects/data/acfhp_trees/for_paml.trees'
directory_out = '/home/alina_grf/progprojects/data/paml_out/'


def run_paml(infile):
    full_path_align = os.path.join(directory_in, infile)
    full_path_out = os.path.join(directory_out, os.path.splitext(infile)[0]+".out")
    print("full_path_out", full_path_out)
    cml = codeml.Codeml(
        alignment=full_path_align,
        tree=tree_file,
        out_file=full_path_out,
        working_dir=directory_out,
    )
    #print(os.path.join(directory_out, 'codeml.ctl'))
    #cml.read_ctl_file(os.path.join(directory_out, 'codeml.ctl'))
    cml.set_options(noisy=9)
    cml.set_options(verbose=1)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=2)
    cml.set_options(CodonFreq=2)
    cml.set_options(model=0)
    cml.set_options(NSsites=[0])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa = 0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=0)
    cml.set_options(omega=.4)
    cml.set_options(fix_alpha=1)
    cml.set_options(alpha=0.)
    cml.set_options(Malpha=0)
    cml.set_options(ncatG=3)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff = .5e-6)
    cml.set_options(cleandata = 0)
    cml.set_options(method=0)
    #cml.print_options()
    try:
        cml.run(command="/home/alina_grf/BIOTOOLS/paml4.9j/bin/codeml", verbose = True)
    except:
        logging.exception("sys.exc_info() {0}, outfile {1}".format(sys.exc_info(), full_path_align))


def parse_dir():
    for infile in os.listdir(directory_in):
        yield infile


if __name__ == '__main__':
    for infile in parse_dir():
        run_paml(infile)

