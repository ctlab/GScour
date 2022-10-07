#!/usr/bin/env python
import logging
import os
import argparse
from Bio import SeqIO
import re

NOT_EQUAL_LENGTH = list()
BROKEN_SPECIES = list()
NOT_MULTIPLE_OF_THREE = list()
EDITED_MULT_OF_THREE = list()
BROKEN_FILES = list()
CORRECT_FILES = list()
BROKEN_FOLDER = "broken_length_f2p"
BROKEN_SPECIES_FOLDER = "broken_species_f2p"
LOG_FILE = "fasta2paml.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
"""
The script consists of two stage:
1. Converting fasta format nucleotide codon sequences (from input directory) to philip-sequential format (to output 
directory)
2. Converting philip-sequential format to specific philip format required by PAML:
In resulting out_dir:  directory "group_id" with folders "file_name" with file_name.phy file for PAML.

For example:
$ cd in_dir
$ ls */
12/:            23/:           12345/:
1.fasta         4055.fasta     2031.fasta 
                3010.fasta     2.fasta

result:
$ cd out_dir
$ ls */
12/:            23/:            12345/:
1/:             4055/:          2031/:   
1.phy           4055.phy        2031.phy
                3010/:          2/:
                3010.phy        2.phy    
"""


def parse_dir_out_gblocks(directory, extens):
    """ parse directory with files out of Gblocks
        .fas-gb can be change just to .fas """
    for personal_folder in os.scandir(directory):
        if os.path.isdir(personal_folder):
            for infile in os.scandir(personal_folder):
                if infile.name.split('.')[-1] == extens:
                    yield personal_folder.name, infile.name


def parse_phylip_dir(in_dir):
    """ parse directory with .phylip files"""
    for personal_folder in os.scandir(in_dir):
        if os.path.isdir(personal_folder):
            for infile in os.scandir(personal_folder):
                if infile.name.split('.')[-1] == 'phylip':
                    yield personal_folder.name, infile.name


def fasta2phylip(personal_folder, infile, folder_in, folder_out):
    """ converting fasta to philip-sequential format"""
    # file_number = re.search(r'(\d+)\.', infile).group(1)
    file_name = infile.split('.')[0]
    infile_path = os.path.join(folder_in, personal_folder, infile)
    outfile_path = os.path.join(folder_out, personal_folder, "{}.{}".format(file_name, "phylip"))
    if not os.path.isdir(os.path.join(folder_out, personal_folder)):
        os.makedirs(os.path.join(folder_out, personal_folder))
    try:
        with open(infile_path, 'r') as input:
            with open(outfile_path, 'w') as output:
                alignments = SeqIO.parse(input, "fasta")
                SeqIO.write(alignments, output, "phylip-sequential")
        logging.info('phylip-sequential format file {} has been recorded'.format(outfile_path))
    except BaseException as err:
        logging.exception("Infile {}, - Unexpected error: {}".format(infile, err))


def chunks(s, n):
    """Produce `n`-character chunks from `s`."""
    for start in range(0, len(s), n):
        if start >= len(s) - n:
            yield s[start:start + n]
        else:
            yield s[start:start + n]+"\n"


def check_lengths(lengths, file_name, species, group):
    """ check: - the lengths of all sequences in one file of the same length
               - number of sequences in one file more or equal then group, less or equal then species
               - length of sequence  a multiple of three
        replace completely broken files to BROKEN_FOLDER
    """
    if all(x == lengths[0] for x in lengths) and group <= len(lengths) <= species and lengths[0] % 3 == 0:
        logging.info("all seq lengths are equal and are the multiple of 3")
        global CORRECT_FILES
        CORRECT_FILES.append(file_name)
    elif not all(x == lengths[0] for x in lengths):
        global NOT_EQUAL_LENGTH
        NOT_EQUAL_LENGTH.append(file_name)
    elif not group <= len(lengths) <= species:
        global BROKEN_SPECIES
        BROKEN_SPECIES.append(file_name)
    elif lengths[0] % 3 != 0:
        global NOT_MULTIPLE_OF_THREE
        NOT_MULTIPLE_OF_THREE.append(file_name)
    else:
        logging.warning("seq lengths are not equal or wrong number of species: {}".format(str(len(lengths))))
        BROKEN_FILES.append(file_name)


def write_target_phy_file(line_edited, target_file):
    for chunk in chunks(line_edited, 60):
        target_file.write(chunk)


def phylip2paml(folder_out, species_folder, source_file_name, species, group):
    """ converting philip-sequential to specific philip format for paml """
    file_name = source_file_name.split('.')[-1]
    target_file_path = os.path.join(folder_out, species_folder, file_name, '{}.{}'.format(file_name, "phy"))
    source_file_path = os.path.join(folder_out, species_folder, source_file_name)
    if not os.path.isdir(os.path.join(folder_out, species_folder, file_name)):
        os.makedirs(os.path.join(folder_out, species_folder, file_name))
    lengths = list()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file_path:
            for line in source_file_path:
                if re.search(r"\d\s(\d+)", line):
                    target_file.write(line)
                elif re.search(r"\d\s{9}", line):
                    """                
                    - insert two spaces and \n instead of 9 after name of sequence
                    - split string on lines by 60 character per line
                    """
                    name_of_seq_9spaces = (re.search(r"(\d)\s{9}", line)).group()
                    name_of_seq = (re.search(r"(\d)", line)).group()
                    target_file.write(name_of_seq + '\n')

                    line_edited = re.sub(name_of_seq_9spaces, "", line)
                    lengths.append(len(line_edited.rstrip()))  # length except \n character
                    write_target_phy_file(line_edited, target_file)

    check_lengths(lengths, file_name, species, group)
    logging.info('changing for paml .phy format file {} has been recorded'.format(target_file_path))


def replace_broken_files(directory_out):
    global BROKEN_FOLDER
    global BROKEN_SPECIES
    if BROKEN_FILES:
        os.makedirs(BROKEN_FOLDER)
        for folder in BROKEN_FILES:
            os.replace(os.path.join(directory_out, folder), os.path.join(BROKEN_FOLDER, folder))
    if BROKEN_SPECIES:  # TODO: testing
        os.makedirs(BROKEN_SPECIES_FOLDER)
        for folder in BROKEN_SPECIES:
            os.replace(os.path.join(directory_out, folder), os.path.join(BROKEN_SPECIES_FOLDER, folder))


def main(folder_in, extension, folder_out, species, group):
    for personal_folder, infile in parse_dir_out_gblocks(folder_in, extension):
        fasta2phylip(personal_folder, infile, folder_in, folder_out)
    for species_folder, phylip_file in parse_phylip_dir(folder_out):
        phylip2paml(folder_out, species_folder, phylip_file, species, group)
        seq_philip_file = os.path.join(folder_out, species_folder, phylip_file)
        os.remove(seq_philip_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with fasta files sorted by separated folders', nargs='?',
                        required=True)
    parser.add_argument('--e', help='Extension of initial FASTA files', nargs='?',
                        default='fas-gb')
    parser.add_argument('--o', help='Path to the folder with result philip files, if it does not exist, it will be'
                                    ' created automatically', nargs='?', required=True)
    parser.add_argument('--species', help='Number of species', nargs='?', required=True)
    parser.add_argument('--group', help='Minimal size of species group', nargs='?', required=True)
    args = parser.parse_args()
    in_dir = args.i
    out_dir = args.o
    ext = args.e
    logging.info("options in use: input folder {}, extension of FASTA {}, output folder {}, number of species {}, "
                 "number of groups {}".format(in_dir, ext, out_dir, args.species, args.group))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    try:
        main(in_dir, ext, out_dir, int(args.species), int(args.group))
        logging.warning("BROKEN_FILES {}:{}".format(len(BROKEN_FILES), BROKEN_FILES))
        if BROKEN_FILES or BROKEN_SPECIES:
            replace_broken_files(out_dir)
        logging.warning("NOT_EQUAL_LENGTH {}:{}".format(len(NOT_EQUAL_LENGTH), NOT_EQUAL_LENGTH))
        logging.warning("BROKEN_SPECIES {}:{}".format(len(BROKEN_SPECIES), BROKEN_SPECIES))
        logging.warning("NOT_MULTIPLE_OF_THREE {}:{}".format(len(NOT_MULTIPLE_OF_THREE), NOT_MULTIPLE_OF_THREE))
        logging.warning("EDITED_MULT_OF_THREE {}:{}".format(len(EDITED_MULT_OF_THREE), EDITED_MULT_OF_THREE))
        logging.info("CORRECT_FILES {}:{}".format(len(CORRECT_FILES), CORRECT_FILES))
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
    logging.info("The work has been completed")
