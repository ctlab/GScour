#!/usr/bin/sudo python
import logging
import os
import argparse
import traceback

from Bio import SeqIO
import re

NOT_EQUAL_LENGTH = list()
NOT_NEEDED_SPECIES = list()
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
In resulting outfolder:  directory "group_id" with folders "file_name" with file_name.phy file for PAML.

For example:
$ cd infolder
$ ls */
12/:            23/:           12345/:
1.fasta         4055.fasta     2031.fasta 
                3010.fasta     2.fasta

result:
$ cd outfolder
$ ls */
12/:            23/:            12345/:
1/:             4055/:          2031/:   
1.phy           4055.phy        2031.phy
                3010/:          2/:
                3010.phy        2.phy    
"""


def parse_dir_out_gblocks(infolder):
    """ parse directory with files out of Gblocks
        .fas-gb can be change just to .fas """
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'fas-gb':
                    yield personal_folder.name, infile


def parse_phylip_dir(infolder):
    """ parse directory with .phylip files"""
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phylip':
                    yield personal_folder.name, infile


def fasta2phylip(personal_folder, infile, folder_in, folder_out):
    """ converting fasta to philip-sequential format"""
    file_number = re.search(r'(\d+)\.', infile).group(1)
    infile_path = os.path.join(folder_in, personal_folder, infile)
    outfile_path = os.path.join(folder_out, personal_folder, "{}.{}".format(file_number, "phylip"))
    if not os.path.isdir(os.path.join(folder_out, personal_folder)):
        os.makedirs(os.path.join(folder_out, personal_folder))
    try:
        with open(infile_path, 'r') as input:
            with open(outfile_path, 'w') as output:
                alignments = SeqIO.parse(input, "fasta")
                SeqIO.write(alignments, output, "phylip-sequential")
        logging.info('phylip-sequential format file {} has been recorded'.format(outfile_path))
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
        global NOT_NEEDED_SPECIES
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


def phylip2paml(folder_out, species_folder, source_file_name, species, group):
    """ converting philip-sequential to specific philip format for paml """
    file_number = re.search(r'(\d+)\.', source_file_name).group(1)
    target_file_path = os.path.join(folder_out, species_folder, file_number, '{}.{}'.format(file_number, "phy"))
    source_file_path = os.path.join(folder_out, species_folder, source_file_name)
    if not os.path.isdir(os.path.join(folder_out, species_folder, file_number)):
        os.makedirs(os.path.join(folder_out, species_folder, file_number))
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
                    lengths.append(len(line_edited.rstrip())) # length except \n character
                    write_target_phy_file(line_edited, target_file)

    check_lengths(lengths, file_number, species, group)
    logging.info('changing for paml and SWAMP .phy format file {} has been recorded'.format(target_file_path))


def replace_broken_files(directory_out):
    global BROKEN_FOLDER
    global NOT_NEEDED_SPECIES
    if BROKEN_FILES:
        os.makedirs(BROKEN_FOLDER)
        for folder in BROKEN_FILES:
            os.replace(os.path.join(directory_out, folder), os.path.join(BROKEN_FOLDER, folder))
    if NOT_NEEDED_SPECIES:  # TODO: testing
        os.makedirs(NOT_NEEDED_SPECIES_FOLDER)
        for folder in NOT_NEEDED_SPECIES:
            os.replace(os.path.join(directory_out, folder), os.path.join(NOT_NEEDED_SPECIES_FOLDER, folder))


def main(folder_in, folder_out, species, group):
    for personal_folder, infile in parse_dir_out_gblocks(folder_in):
        fasta2phylip(personal_folder, infile, folder_in, folder_out)
    for species_folder, phylip_file in parse_phylip_dir(folder_out):
        phylip2paml(folder_out, species_folder, phylip_file, species, group)
        seq_philip_file = os.path.join(folder_out, species_folder, phylip_file)
        os.remove(seq_philip_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with fasta files sorted by separated folders', nargs='?')
    parser.add_argument('--outfolder', help='Path to the folder with result philip files', nargs='?')
    parser.add_argument('--species', help='Number of species', nargs='?')
    parser.add_argument('--group', help='Minimal size of species group', nargs='?')
    args = parser.parse_args()
    outfolder = args.outfolder
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    try:
        main(args.infolder, outfolder, int(args.species), int(args.group))
        logging.warning("BROKEN_FILES {}:{}".format(len(BROKEN_FILES), BROKEN_FILES))
        if BROKEN_FILES or NOT_NEEDED_SPECIES:
            replace_broken_files(outfolder)
        logging.warning("NOT_EQUAL_LENGTH {}:{}".format(len(NOT_EQUAL_LENGTH), NOT_EQUAL_LENGTH))
        logging.warning("NOT_NEEDED_SPECIES {}:{}".format(len(NOT_NEEDED_SPECIES), NOT_NEEDED_SPECIES))
        logging.warning("NOT_MULTIPLE_OF_THREE {}:{}".format(len(NOT_MULTIPLE_OF_THREE), NOT_MULTIPLE_OF_THREE))
        logging.warning("EDITED_MULT_OF_THREE {}:{}".format(len(EDITED_MULT_OF_THREE), EDITED_MULT_OF_THREE))
    except BaseException as e:
        logging.info("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))
    logging.info("The work has been completed")
