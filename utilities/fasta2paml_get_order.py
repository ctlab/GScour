#!/usr/bin/env python
import logging
import os
import argparse
from itertools import *
import subprocess
from Bio.Phylo.PAML import codeml
from Bio import SeqIO
import re
import shutil

NOT_EQUAL_LENGTH = list()
BROKEN_SPECIES = list()
NOT_MULTIPLE_OF_THREE = list()
EDITED_MULT_OF_THREE = list()
BROKEN_FILES = list()
LOG_FILE = "fasta2paml_get_order.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
"""
This script trying to get right order of sequences for paml input file
------------------------------------------------
The script consists of two stage:
1. Converting fasta format nucleotide codon sequences (from input directory) to philip-sequential format (to output 
directory) in accordance with specific order 
(which is searched in the file in_dir/species_folder/species_folder_name.order) required for the paml
2. Converting philip-sequential format to specific philip format required by PAML:
In resulting out_dir:  directory "group_id" with folders "file_name" with file_name.phy file for PAML.

For example:
$ cd in_dir
$ ls */
12/:            23/:           12345/:
1.fasta         4055.fasta     2031.fasta 
                3010.fasta     2.fasta
$ cd trees/
$ ls */
12.tree
23.tree
12345.tree
                
$ cd order_dir
$ ls */
12.order    23.order    12345.order

result:
$ cd out_dir
$ ls */
12/:            23/:            12345/:
1/:             4055/:          2031/:   
1.phy           4055.phy        2031.phy
                3010/:          2/:
                3010.phy        2.phy    
"""


def get_order(folder_in, species_folder):
    """ parse directory with .order files
    which contain specific order of seqs appearance
    if there is no order file just return name of species folder
    separated by comma as order_string"""
    for root, dirs, files in os.walk(folder_in):
        order_file_name = '{}.{}'.format(species_folder.name, 'order')
        if order_file_name in files:
            with open(os.path.join(folder_in, order_file_name), 'r') as f:
                order_string = f.read()
                return order_string.rstrip()
        else:
            ll = list(species_folder.name)
            order_string = ','.join(ll)
            return order_string


def get_infile_and_order(folder_in, folder_order, logger):
    """ parse directory with files out of Gblocks
        'fas-gb' can be change just to .fas"""
    order_string = get_order(folder_order, folder_in)  # e.g 1,2,3,5,6
    print("order string", order_string)
    if not order_string:
        logger.warning("Please check .order file for {}".format(folder_in.name))
        return
    for infile in os.scandir(folder_in):
        if infile.name.split('.')[-1] == 'fas-gb':
            print('return fas-gb', infile.name)
            yield infile.name, order_string.split(',')


def get_phylip_file(in_folder, folder_name):
    """ parse directory with .phylip files"""
    for species_folder in os.scandir(in_folder):
        if os.path.isdir(species_folder) and species_folder.name == folder_name:
            for infile in os.listdir(species_folder):
                if infile.split('.')[-1] == 'phylip':
                    print("return phylip", infile)
                    return infile


def fasta2phylip(personal_folder, infile, orders, folder_in, folder_out, logger):
    """ converting fasta to philip-sequential format"""
    logger.info("start converting fasta2phylip")
    # file_number = re.search(r'(\d+)\.', infile).group(1) # TODO: make it possible to name files with letters
    file_name = infile.split('.')[0]
    infile_path = os.path.join(folder_in, personal_folder, infile)
    outfile_path = os.path.join(folder_out, personal_folder, "{}.{}".format(file_name, "phylip"))
    print("outfile_path", outfile_path)
    if not os.path.isdir(os.path.join(folder_out, personal_folder)):
        os.makedirs(os.path.join(folder_out, personal_folder))
    try:
        with open(infile_path, 'r') as input_file:
            with open(outfile_path, 'a') as output_file:
                alignments = SeqIO.parse(input_file, "fasta")
                ordering_alignments = list()
                alignments = list(alignments)
                for name_of_seq in orders:
                    for align in alignments:
                        if align.name == name_of_seq:
                            ordering_alignments.append(align)
                            break
                SeqIO.write(ordering_alignments, output_file, "phylip-sequential")
                return outfile_path
    except BaseException as e:
        logger.error("BaseException: {} for file {}".format(e, infile))


def chunks(s, n):
    """Produce `n`-character chunks from `s`."""
    for start in range(0, len(s), n):
        # yield s[start:start + n]+'\n'
        if start >= len(s) - n:
            yield s[start:start + n]
        else:
            yield s[start:start + n]  # +"\n"


def check_lengths(lengths, species_folder, file_name, species, group, logger):
    """ check: - the lengths of all sequences in one file are of the same length
               - number of sequences in one file more or equal then group number, less or equal then species number
               - length of sequence is multiple of three
        replace completely broken files to BROKEN_FOLDER
    """
    if all(x == lengths[0] for x in lengths) and group <= len(lengths) <= species and lengths[0] % 3 == 0:
        logging.info("all seq lengths are equal and are the multiple of 3")
        global CORRECT_FILES
        CORRECT_FILES.append(file_name)
    elif not all(x == lengths[0] for x in lengths):
        global NOT_EQUAL_LENGTH
        NOT_EQUAL_LENGTH.append('{}/{}'.format(species_folder, file_name))
    elif not group <= len(lengths) <= species:
        global BROKEN_SPECIES
        BROKEN_SPECIES.append('{}/{}'.format(species_folder, file_name))
    elif lengths[0] % 3 != 0:
        global NOT_MULTIPLE_OF_THREE
        NOT_MULTIPLE_OF_THREE.append('{}/{}'.format(species_folder, file_name))
    else:
        logger.warning("seq lengths are not equal or wrong number of species: {}".format(str(len(lengths))))
        BROKEN_FILES.append('{}/{}'.format(species_folder, file_name))


def write_target_phy_file(line_edited, target_file):
    for chunk in chunks(line_edited, 60):
        target_file.write(chunk)


def phylip2paml(folder_out, species_folder, source_file_name, species, group, logger):
    """ converting philip-sequential to specific philip format for paml """
    logging.info("start converting phylip2paml")
    print("source_file_path:", folder_out, species_folder, source_file_name)
    source_file_path = source_file_name # os.path.join(folder_out, species_folder, source_file_name)
    file_name = (source_file_path.split('/')[-1]).split('.')[0]
    # file_number = re.search(r'(\d+)\.', source_file_name).group(1)
    target_file_path = os.path.join(folder_out, species_folder, file_name, '{}.{}'.format(file_name, "phy"))
    print("target_file_path", target_file_path)
    if not os.path.isdir(os.path.join(folder_out, species_folder, file_name)):
        os.makedirs(os.path.join(folder_out, species_folder, file_name))
    lengths = list()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file_path:
            for line in source_file_path:
                if re.search(r"\d\s(\d+)", line):
                    target_file.write(line)
                elif re.search(r"\d+\s{9}", line):
                    """                
                    - insert two spaces and \n instead of 9 after name of sequence
                    - split string on lines by 60 character per line
                    """
                    name_of_seq_9spaces = (re.search(r"(\d+)\s{9}", line)).group()
                    name_of_seq = (re.search(r"(\d+)", line)).group()
                    target_file.write(name_of_seq + '\n')

                    line_edited = re.sub(name_of_seq_9spaces, "", line)
                    # lengths.append(len(line_edited.rstrip()))  # length except \n character
                    lengths.append(len(line_edited[:-1]))  # length except \n character
                    write_target_phy_file(line_edited, target_file)

    check_lengths(lengths, species_folder, file_name, species, group, logger)
    logger.info('changing for paml and SWAMP .phy format file {} has been recorded'.format(target_file_path))
    return target_file_path


def replace_broken_files(directory_out):
    broken_length_folder = os.path.join(directory_out, "broken_length_files")
    not_needed_species_folder = os.path.join(directory_out, "not_needed_species")
    if BROKEN_FILES:
        for folder in BROKEN_FILES:
            shutil.move(os.path.join(directory_out, folder), os.path.join(broken_length_folder, folder))
    if BROKEN_SPECIES:
        for folder in BROKEN_SPECIES:
            shutil.move(os.path.join(directory_out, folder), os.path.join(not_needed_species_folder, folder))
            # os.remove(os.path.join(directory_out, folder)) # TODO: shutil not delete source, just leave empty


def get_tree_path(trees_folder, species_folder_name):
    for tree in os.scandir(trees_folder):
        if tree.name.split('.')[0] == species_folder_name:
            return tree.name


def get_input_items(folder_in, trees_folder, folder_name, logger):
    """ parse root folder with files for paml
    parse tree_folder to get appropriate tree
    :param folder_name:
    :param trees_folder:
    :param folder_in: 
    :param logger: """
    for species_folder in os.scandir(folder_in):
        print("folder_in", folder_in, "species_folder", species_folder.name,"folder_name", folder_name)
        if os.path.isdir(species_folder) and species_folder.name == folder_name:
            tree_name = get_tree_path(trees_folder, species_folder.name)
            print("tree_name", tree_name)
            if not tree_name:
                logger.warning("There is no tree file for group {}".format(species_folder.name))
                return
            tree_path = os.path.join(trees_folder, tree_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.scandir(item):
                        if infile.name.split('.')[-1] == 'phy': # and infile.split('.')[0].isnumeric():
                            return item.name, infile.name, tree_path


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
    cml.set_options(seqtype=1)  # 1:codons; 2:AAs; 3:codons-->AAs
    cml.set_options(CodonFreq=2)  # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    cml.set_options(clock=0)  # 0:no clock, 1:global clock; 2:local clock; 3:TipDate
    cml.set_options(aaDist=0)
    cml.set_options(model=0)  # models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches
    cml.set_options(NSsites=[0])  # * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
    # 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
    # 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
    # 13:3normal>0
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)  # 1: kappa fixed, 0: kappa to be estimated
    cml.set_options(kappa=2)  # initial or fixed kappa
    cml.set_options(fix_omega=0)  # 1:omega fixed, 0:omega to be estimated
    cml.set_options(omega=1)  # initial or fixed omega, for codons or codon-based AAs
    cml.set_options(fix_alpha=1)  # 0: estimate gamma shape parameter; 1: fix it at alpha
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


def run_codeml(folder_in, species_folder, item_folder, infile, phylogeny_tree_path, exec_path, time_out, logger):
    item_folder_path = os.path.join(folder_in, species_folder, item_folder)
    file_name = infile.split('.')[0]
    # try:
    #     file_number = (re.search(r"(\d+).phy", infile)).group(1)
    # except AttributeError as err:
    #     logger.info("Please check file .phy {}: cause an error {}".format(item_folder_path, err.args))
    #     return
    os.chdir(item_folder_path)
    logger.info("Codeml: working with {}".format(file_name))
    infile_path = os.path.join(item_folder_path, infile)
    set_one_ratio_model(infile_path, phylogeny_tree_path, item_folder_path)
    try:
        p = subprocess.Popen(exec_path, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        # return_code = p.wait(timeout=time_out)  # Timeout in seconds
        stdout, stderr = p.communicate(input=b'\n', timeout=time_out)
        # p.communicate(input='\n')
        logger.info("OK paml for file number {}, stderr: {}".format(file_name, stderr))
        # infile_path_copy = os.path.join(item_folder_path, "{}_{}.{}".format(infile.split('.')[0], 'copy',
        #                                                                     infile.split('.')[1]))
        # shutil.copyfile(infile_path, infile_path_copy)
        # if return_code:
        #     return False
        # return True  # if return_code == 0
    except subprocess.SubprocessError as err:
        logger.info("SubprocessError: {}, \nerr.args".format(err, err.args))
    except subprocess.CalledProcessError as err:
        logger.info("CalledProcessError: {}, \nerr.stderr: {}\nerr.args\n: {}".format(err, err.stderr, err.args))
            # exit(1)
    except subprocess.TimeoutExpired as err:
        p.kill()
        file_id = "{}/{}".format(species_folder, item_folder)
        logger.warning("Killed {}, {}, try to increase timeout or launch codeml manually and check,"
                       .format(file_id, err.args))
    else:
        return True


def wright_order_file(folder_in, species_folder, order_list, logger):
    order_file_path = os.path.join(folder_in, '{}.{}'.format(species_folder, 'order'))
    order_string = ','.join(order_list)
    if os.path.isfile(order_file_path):
        open(order_file_path, 'w').close()  # clean file
    with open(order_file_path, 'w') as f:
        f.write(order_string)
        logger.info("file with right order {} recorded".format(order_file_path))


def main(folder_in, folder_order, folder_trees, exec_path, time_out, folder_out, species, group, logger):
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder): # and species_folder.name.isdigit():
            guess = False
            for infile_fasta, order_list in get_infile_and_order(species_folder, folder_order, logger):
                if guess:
                    break
                print("order_list", order_list)
                # list_of_orders = list(permutations(order_list))
                # print("list_of_orders", list_of_orders)
                while not guess:
                    for order in permutations(order_list):
                        phylip_file = fasta2phylip(species_folder, infile_fasta, order, folder_in, folder_out, logger)
                        # phylip_file = get_phylip_file(folder_out, species_folder)
                        infile_phy = phylip2paml(folder_out, species_folder.name, phylip_file, species, group, logger)
                        seq_philip_file = os.path.join(folder_out, species_folder.name, phylip_file)
                        os.remove(seq_philip_file)
                        logger.info("remove {}".format(seq_philip_file))
                        if not get_input_items(folder_out, folder_trees, species_folder.name, logger):
                            break
                        item_folder, infile_phy, tree_path = get_input_items(folder_out, folder_trees, species_folder.name,
                                                                             logger)
                        logger.info("run codeml")
                        if not run_codeml(folder_out, species_folder.name, item_folder, infile_phy, tree_path, exec_path,
                                          time_out,
                                          logger):
                            logger.info('FAIL!Wrong order {} for species folder {}'.format(order, species_folder.name))
                            continue
                            # try:
                            #     list_of_orders.pop(0)
                            #     order_list = list_of_orders[0]
                            # except IndexError:
                            #     logger.warning("list of orders is {} empty, try to check codeml manually".format(list_of_orders))
                            #     guess = True
                        else:
                            guess = True
                            logger.info('GUESS!Right order for species folder {}={}'.format(species_folder.name,
                                                                                            order_list))
                            wright_order_file(folder_order, species_folder.name, order_list, logger)
                            input_dir = os.path.join(folder_in, species_folder)
                            shutil.rmtree(input_dir)
                            logging.info("Input dir {} with fasta file deleted".format(input_dir))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the input folder with fasta files sorted by separated folders', nargs='?',
                        required=True)
    parser.add_argument('--order', help='Path to the folder with .order files for each folder in input, this folder '
                                        'can be initially empty', nargs='?',
                        required=True)
    parser.add_argument('--tree', help='Path to the folder with trees for paml, file of tree must be named as species '
                                       'folder', nargs='?', required=True)
    parser.add_argument('--e', help='Path to the codeml executable', nargs='?', default='codeml')
    parser.add_argument('--timeout', help='Timeout for codeml in seconds, default=120', nargs='?', default='120')
    parser.add_argument('--o', help='Path to the folder with result philip files, can be created automatically',
                        nargs='?', required=True)
    parser.add_argument('--species', help='Number of species', nargs='?', required=True)
    parser.add_argument('--group', help='Minimal size of species group', nargs='?', required=True)
    args = parser.parse_args()
    out_dir = args.o
    executable_path = args.e
    timeout = int(args.timeout)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    child_logger = logging.getLogger(__name__)
    child_logger.addHandler(logging.FileHandler(LOG_FILE))
    child_logger.setLevel(10)
    try:
        main(args.i, args.order, args.tree, executable_path, timeout, out_dir, int(args.species), int(args.group),
             child_logger)
        logging.warning("BROKEN_FILES {}:{}".format(len(BROKEN_FILES), BROKEN_FILES))
        if BROKEN_FILES or BROKEN_SPECIES:
            replace_broken_files(out_dir)
        child_logger.warning("NOT_EQUAL_LENGTH {}:{}".format(len(NOT_EQUAL_LENGTH), NOT_EQUAL_LENGTH))
        child_logger.warning("BROKEN_SPECIES number {}:{}".format(len(BROKEN_SPECIES), BROKEN_SPECIES))
        child_logger.warning("NOT_MULTIPLE_OF_THREE {}:{}".format(len(NOT_MULTIPLE_OF_THREE), NOT_MULTIPLE_OF_THREE))
        child_logger.warning("EDITED_MULT_OF_THREE {}:{}".format(len(EDITED_MULT_OF_THREE), EDITED_MULT_OF_THREE))
        child_logger.info("CORRECT_FILES {}:{}".format(len(CORRECT_FILES), CORRECT_FILES))
    except BaseException as e:
        child_logger.exception("Unexpected error: {}".format(e))
        raise e
    child_logger.info("The work has been completed")
