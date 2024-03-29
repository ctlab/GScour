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
BROKEN_FOLDER = "broken_length_files(fasta2paml)"
NOT_NEEDED_SPECIES_FOLDER = "not_needed_species(f2p)"
LOG_FILE = "fasta2paml.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
"""
The script consists of two stage:
1. Converting fasta format nucleotide codon sequences (from infolder) to philip-sequential format (to outfolder)
2. Converting philip-sequential format to specific philip format required by PAML:
In resulting outfolder: number_of_file.philip and directory "number_of_file" with .phy file for PAML.

It is the first version of file where no division into groups was implied. But there some check:
group <= len(number of sequences in one file) <= species

Script works, but results which are not divided into separate folders by groups
would be hard for further processing (with paml, swamp)
"""


def parse_dir_out_gblocks(infolder):
    """ parse directory with files out of Gblocks
        'fas-gb' can be change just to .fas extension """
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fas-gb':
            yield os.path.join(infolder, infile)


def parse_phylip_dir(infolder):
    """ parse directory with .phylip files"""
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'phylip':
            yield os.path.join(infolder, infile)


def fasta2phylip(infile, folder_out):
    """ converting fasta to philip-sequential format"""
    outfile = os.path.join(folder_out, "{}.{}".format(re.search(r'\/(\d+)\.', infile).group(1), "phylip"))
    try:
        with open(infile, 'r') as input:
            with open(outfile, 'w') as output:
                alignments = SeqIO.parse(input, "fasta")
                SeqIO.write(alignments, output, "phylip-sequential")
        logging.info('phylip-sequential format file {} has been recorded'.format(outfile))
    except BaseException as e:
        logging.error("BaseException: {} for file {}".format(e, infile))


def chunks(s, n):
    """Produce `n`-character chunks from `s`."""
    for start in range(0, len(s), n):
        if start >= len(s) - n:
            yield s[start:start + n]
        else:
            yield s[start:start + n]+"\n"


def check_lengths(lengths, file_number, species, group):
    """ check: - the lengths of all sequences in one file of the same length
               - number of sequences in one file more or equal then group, less or equal then species
               - length of sequence  a multiple of three
        replace completely broken files to BROKEN_FOLDER
    """
    if all(x == lengths[0] for x in lengths) and group <= len(lengths) <= species and lengths[0] % 3 == 0:
        logging.info("all seq lengths are equal, confirm of number of species")
    elif not all(x == lengths[0] for x in lengths):
        global NOT_EQUAL_LENGTH
        NOT_EQUAL_LENGTH.append(file_number)
    elif not group <= len(lengths) <= species:
        global BROKEN_SPECIES
        NOT_NEEDED_SPECIES.append(file_number)
    elif lengths[0] % 3 != 0:
        global NOT_MULTIPLE_OF_THREE
        NOT_MULTIPLE_OF_THREE.append(file_number)
    else:
        logging.warning("seq lengths are not equal or wrong number of species: {}".format(str(len(lengths))))
        BROKEN_FILES.append(file_number)


def write_target_phy_file(line_edited, target_file):
    for chunk in chunks(line_edited, 60):
        target_file.write(chunk)


def phylip2paml(source_file_path, species, group):
    """ converting philip-sequential to specific philip format for paml """
    file_number = re.search(r'(\d+)\.', source_file_path).group(1)
    personal_folder = os.path.join(os.path.split(source_file_path)[0], '{}'.format(file_number))
    target_file_path = os.path.join(personal_folder, '{}.{}'.format(file_number, "phy"))
    if not os.path.isdir(personal_folder):
        os.makedirs(personal_folder)
    lengths = list()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
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

    check_lengths(lengths, file_number, species, group)
    logging.info('changing for paml and SWAMP .phy format file {} has been recorded'.format(target_file_path))


def replace_broken_files(directory_out):
    global BROKEN_FOLDER
    global BROKEN_SPECIES
    if BROKEN_FILES:
        os.makedirs(BROKEN_FOLDER)
        for folder in BROKEN_FILES:
            os.replace(os.path.join(directory_out, folder), os.path.join(BROKEN_FOLDER, folder))
    if NOT_NEEDED_SPECIES:  # TODO: testing
        os.makedirs(NOT_NEEDED_SPECIES_FOLDER)
        for folder in NOT_NEEDED_SPECIES:
            os.replace(os.path.join(directory_out, folder), os.path.join(NOT_NEEDED_SPECIES_FOLDER, folder))


def main(folder_in, folder_out, species, group):
    for infile in parse_dir_out_gblocks(folder_in):
        fasta2phylip(infile, folder_out)
    for phylip_file in parse_phylip_dir(folder_out):
        phylip2paml(phylip_file, species, group)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with fasta files', nargs='?', required=True)
    parser.add_argument('--o', help='Path to the folder with result philip files', nargs='?', required=True)
    parser.add_argument('--species', help='Number of species', nargs='?', required=True)
    parser.add_argument('--group', help='Minimal size of species group', nargs='?', required=True)
    args = parser.parse_args()
    out_folder = args.o
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    try:
        main(args.i, out_folder, int(args.species), int(args.group))
        logging.warning("BROKEN_FILES {}:{}".format(len(BROKEN_FILES), BROKEN_FILES))
        if BROKEN_FILES or BROKEN_SPECIES:
            replace_broken_files(out_folder)
        logging.warning("NOT_EQUAL_LENGTH {}:{}".format(len(NOT_EQUAL_LENGTH), NOT_EQUAL_LENGTH))
        logging.warning("NOT_NEEDED_SPECIES {}:{}".format(len(BROKEN_SPECIES), BROKEN_SPECIES))
        logging.warning("NOT_MULTIPLE_OF_THREE {}:{}".format(len(NOT_MULTIPLE_OF_THREE), NOT_MULTIPLE_OF_THREE))
        logging.warning("EDITED_MULT_OF_THREE {}:{}".format(len(EDITED_MULT_OF_THREE), EDITED_MULT_OF_THREE))
    except BaseException as err:
        logging.info("Unexpected error: {}".format(err))
    logging.info("The work has been completed")
