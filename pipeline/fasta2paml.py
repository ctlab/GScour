#!/usr/bin/env python
import logging
import os
import argparse
from Bio import SeqIO
import re

BROKEN_FILES = list()
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
        if start >= len(s) - n:
            yield s[start:start + n]
        else:
            yield s[start:start + n]+"\n"


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
                    line_edited_end = re.sub(r"T\s?G\s?A$", "", line)
                    line_edited_end = re.sub(r"T\s?A\s?G$", "", line)
                    line_edited_end = re.sub(r"T\s?A\s?A$", "", line)
                    repl_9_spaces = (re.search(r"(\d)\s{9}", line)).group()
                    repl_2_spaces = (re.search(r"(\d)", line)).group()
                    target_file.write(repl_2_spaces + '\n')
                    line_edited = re.sub(repl_9_spaces, "", line_edited_end)
                    for chunk in chunks(line_edited, 60):
                        target_file.write(chunk)
                if re.search(r"\d\s(\d+)", line):
                    """
                    decreasing the number of characters because of removing stop codon
                    """
                    repl = (re.search(r"\d\s(\d+)", line)).group(1)
                    repl = str(int(repl) - 3)
                    line_edited = re.sub(r"\d+$", repl+"", line)
                    target_file.write(line_edited)
    if all(x == lengths[0] for x in lengths) and len(lengths) == species:
        logging.info("all seq lengths are equal, confirm of number of species")
    else:
        logging.warning("seq lengths are not equal or wrong number of species: {}".format(str(len(lengths))))
        BROKEN_FILES.append(file_number)
    logging.info('changing for paml and SWAMP .phy format file {} has been recorded'.format(target_file_path))


def main(infolder, outfolder, species):
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
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
    try:
        main(args.infolder, args.outfolder, int(args.species))
        if BROKEN_FILES:
            logging.warning("BROKEN_FILES:{}".format(BROKEN_FILES))
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
