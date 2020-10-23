#!/usr/bin/env python
import logging
import os
import argparse
from Bio import SeqIO
import re

NOT_EQUAL_LENGTH = list()
NOT_NEEDED_SPECIES = list()
NOT_MULTIPLE_OF_THREE = list()
BROKEN_FILES = list()
STOP_CODON_BROKEN_FILES = list()

LOG_FILE = "fasta2paml.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def parse_dir_out_gblocks(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fas-gb':
            yield os.path.join(infolder, infile)


def parse_phylip_dir(infolder):
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'phylip':
            yield os.path.join(infolder, infile)


def fasta2phylip(infile, outfolder):
    outfile = os.path.join(outfolder, "{}.{}".format(re.search(r'\/(\d+)\.', infile).group(1), "phylip"))
    with open(infile, 'r') as input:
        with open(outfile, 'w') as output:
            alignments = SeqIO.parse(input, "fasta")
            SeqIO.write(alignments, output, "phylip-sequential")
    logging.info('phylip-sequential format file {} has been recorded'.format(outfile))


def chunks(s, n):
    """Produce `n`-character chunks from `s`."""
    for start in range(0, len(s), n):
        yield s[start:start + n]
        """
        if start >= len(s) - n:
            yield s[start:start + n]
        else:
            yield s[start:start + n]+"\n"
        """


def phylip2paml(source_file_path, species):
    file_number = re.search(r'(\d+)\.', source_file_path).group(1)
    personal_folder = os.path.join(os.path.split(source_file_path)[0], '{}'.format(file_number))
    target_file_path = os.path.join(personal_folder, '{}.{}'.format(file_number, "phy"))
    if not os.path.isdir(personal_folder):
        os.makedirs(personal_folder)
    lengths = list()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
                if re.search(r"\d\s{9}", line):
                    """
                    - remove stop codon
                    - insert two spaces and \n instead of 9 after name of sequence
                    - split string on lines by 60 character per line
                    """
                    lengths.append(len(line))
                    if re.search(r"T\s?G\s?A$", line):
                        line_edited_end = re.sub(r"T\s?G\s?A$", "", line)
                    if re.search(r"T\s?A\s?G$", line):
                        line_edited_end = re.sub(r"T\s?A\s?G$", "", line)
                    if re.search(r"T\s?A\s?A$", line):
                        line_edited_end = re.sub(r"T\s?A\s?A$", "", line)
                    repl_9_spaces = (re.search(r"(\d)\s{9}", line)).group()
                    repl_2_spaces = (re.search(r"(\d)", line)).group()
                    target_file.write(repl_2_spaces + '\n')
                    try:
                        line_edited = re.sub(repl_9_spaces, "", line_edited_end)
                        for chunk in chunks(line_edited, 60):
                            target_file.write(chunk)
                    except UnboundLocalError:
                        logging.warning("no stop codon in file {}".format(source_file))
                        if file_number not in STOP_CODON_BROKEN_FILES:
                            STOP_CODON_BROKEN_FILES.append(file_number)
                if re.search(r"\d\s(\d+)", line):
                    """
                    decreasing the number of characters because of removing stop codon
                    """
                    repl = (re.search(r"\d\s(\d+)", line)).group(1)
                    repl = str(int(repl) - 3)
                    line_edited = re.sub(r"\d+$", repl+"", line)
                    target_file.write(line_edited)
    if all(x == lengths[0] for x in lengths) and len(lengths) == species and lengths[0] % 3 == 0:
        logging.info("all seq lengths are equal, confirm of number of species")
    elif not all(x == lengths[0] for x in lengths):
        global NOT_EQUAL_LENGTH
        NOT_EQUAL_LENGTH.append(file_number)
    elif not len(lengths) == species:
        global NOT_NEEDED_SPECIES
        NOT_NEEDED_SPECIES.append(file_number)
    elif lengths[0] % 3 != 0:
        global NOT_MULTIPLE_OF_THREE
        NOT_MULTIPLE_OF_THREE.append(file_number)
    else:
        logging.warning("seq lengths are not equal or wrong number of species: {}".format(str(len(lengths))))
        BROKEN_FILES.append(file_number)
    logging.info('changing for paml and SWAMP .phy format file {} has been recorded'.format(target_file_path))


def replace_broken_files(directory_out):
    broken_folder = "broken_length_files(fasta2paml)"
    os.makedirs(broken_folder)
    for folder in BROKEN_FILES:
        os.replace(os.path.join(directory_out, folder), os.path.join(broken_folder, folder))


def main(infolder, outfolder, species):
    for infile in parse_dir_out_gblocks(infolder):
        fasta2phylip(infile, outfolder)
    for phylip_file in parse_phylip_dir(outfolder):
        phylip2paml(phylip_file, species)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--outfolder', help='Path to the folder with output files of prank', nargs='?')
    parser.add_argument('--species', help='Number of species', nargs='?')
    args = parser.parse_args()
    outfolder = args.outfolder
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    try:
        main(args.infolder, outfolder, int(args.species))
        if BROKEN_FILES:
            logging.warning("BROKEN_FILES {}:{}".format(len(BROKEN_FILES), BROKEN_FILES))
            replace_broken_files(outfolder)
        if STOP_CODON_BROKEN_FILES:
            logging.warning("STOP_CODON_BROKEN_FILES:{}".format(STOP_CODON_BROKEN_FILES))
        logging.warning("NOT_EQUAL_LENGTH {}:{}".format(len(NOT_EQUAL_LENGTH), NOT_EQUAL_LENGTH))
        logging.warning("NOT_NEEDED_SPECIES {}:{}".format(len(NOT_NEEDED_SPECIES), NOT_NEEDED_SPECIES))
        logging.warning("NOT_MULTIPLE_OF_THREE {}:{}".format(len(NOT_MULTIPLE_OF_THREE), NOT_MULTIPLE_OF_THREE))
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
