#!/usr/bin/env python
import logging
import os
import argparse
from Bio import AlignIO
import sys
import re

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


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
            alignments = AlignIO.parse(input, "fasta")
            AlignIO.write(alignments, output, "phylip-sequential")
    logging.info('phylip format file {} has been recorded'.format(outfile))


def phylip2paml(source_file_path):
    corrected_file_path = os.path.join(os.path.split(source_file_path)[0], '{}.{}'.format(re.search(r'\/(\d+)\.',
                                                                                 source_file_path).group(1), "phy"))
    with open(corrected_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
                if re.search(r"\d\s{9}", line):
                    line_edited_end = re.sub(r"T\s?G\s?A$", "\n", line)
                    repl_9_spaces = (re.search(r"(\d)\s{9}", line)).group()
                    repl_2_spaces = (re.search(r"(\d)\s{2}", line)).group()
                    line_edited = re.sub(repl_9_spaces, repl_2_spaces, line_edited_end)
                    target_file.write(line_edited)
                else:
                    repl = (re.search(r"\d\s(\d+)", line)).group(1)
                    repl = str(int(repl) - 3)
                    line_edited = re.sub(r"\d+$", repl+"\n", line)
                    target_file.write(line_edited)
    logging.info('phy format file {} has been recorded'.format(corrected_file_path))


def main(infolder, outfolder):
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    for infile in parse_dir_out_gblocks(infolder):
        fasta2phylip(infile, outfolder)
    for phylip_file in parse_phylip_dir(outfolder):
        phylip2paml(phylip_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for prank', nargs='?')
    parser.add_argument('--outfolder', help='Path to the folder with output files of prank', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder, args.outfolder)
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
