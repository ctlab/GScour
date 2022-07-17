#!/usr/bin/env python
import logging
import os
import argparse
import sys
import traceback

import pandas as pd
from Bio import SeqIO
import re
import shutil

NOT_EQUAL_LENGTH = list()
BROKEN_SPECIES = dict()
NOT_MULTIPLE_OF_THREE = list()
BROKEN_LENGTH_FILES = dict()
LOG_FILE = "fasta2paml_ordering.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
"""
The script consists of two stage:
1. Converting fasta format nucleotide codon sequences (from input directory) to philip-sequential format (to output 
directory)
in accordance with specific order (which is searched in the file in_dir/species_folder/species_folder_name.order) 
required for the paml
2. Converting philip-sequential format to specific philip format required by PAML:
In resulting out_dir:  directory "group_id" with folders "file_name" with file_name.phy file for PAML.

For example:
$ cd in_dir
$ ls */
12/:            23/:           12345/:
1.fasta         4055.fasta     2031.fasta 
                3010.fasta     2.fasta
                               
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
    which contain specific order of seqs appearance """
    for infile in os.scandir(folder_in):
        if infile.name == '{}.{}'.format(species_folder, 'order'):
            with open(os.path.join(folder_in, infile.name), 'r') as f:
                order_string = f.read()
                return order_string.rstrip()


def get_infile_and_order(folder_in, folder_order):
    """ parse directory with files out of Gblocks
        'fas-gb' can be change just to .fas"""
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder) and species_folder.name.isdigit():
            order_string = get_order(folder_order, species_folder.name)
            if not order_string:
                logging.warning("Please check .order file for {}{}".format(folder_in, species_folder.name))
                return
            for infile in os.scandir(species_folder):
                if infile.name.split('.')[-1] == 'fas-gb':
                    yield species_folder.name, infile.name, order_string


def parse_phylip_dir(input_folder):
    """ parse directory with .phylip files"""
    for personal_folder in os.scandir(input_folder):
        if os.path.isdir(personal_folder):
            for infile in os.scandir(personal_folder):
                if infile.name.split('.')[-1] == 'phylip':
                    yield personal_folder.name, infile.name


def fasta2phylip(personal_folder, infile_name, order_string, folder_in, folder_out):
    """ converting fasta to philip-sequential format"""
    # file_number = re.search(r'(\d+)\.', infile).group(1)
    file_gene_name = infile_name.split('.')[0]
    infile_path = os.path.join(folder_in, personal_folder, infile_name)
    outfile_path = os.path.join(folder_out, personal_folder, "{}.{}".format(file_gene_name, "phylip"))
    if not os.path.isdir(os.path.join(folder_out, personal_folder)):
        os.makedirs(os.path.join(folder_out, personal_folder))
    try:
        with open(infile_path, 'r') as input_file:
            with open(outfile_path, 'a') as output_file:
                alignments = SeqIO.parse(input_file, "fasta")
                ordering_alignments = list()
                alignments = list(alignments)
                for name_of_seq in order_string.split(','):
                    for align in alignments:
                        if align.name == name_of_seq:
                            ordering_alignments.append(align)
                            break
                SeqIO.write(ordering_alignments, output_file, "phylip-sequential")
        logging.info('phylip-sequential format file {} has been recorded'.format(outfile_path))
    except BaseException as e:
        logging.exception("Infile {}, - Unexpected error: {}".format(infile_name, e))


def chunks(s, n):
    """Produce `n`-character chunks from `s`."""
    for start in range(0, len(s), n):
        # yield s[start:start + n]+'\n'
        if start >= len(s) - n:
            yield s[start:start + n]
        else:
            yield s[start:start + n]  # +"\n"


def check_lengths(lengths, species_folder, file_name, species, group):
    """ check: - the lengths of all sequences in one file are of the same length
               - number of sequences in one file more or equal then group number, less or equal then species number
               - length of sequence is multiple of three
        replace completely broken files to BROKEN_FOLDER
    """
    if all(x == lengths[0] for x in lengths) and group <= len(lengths) <= species and lengths[0] % 3 == 0:
        logging.info("Check length OK: all seq lengths are equal, multiple of 3")
    elif not all(x == lengths[0] for x in lengths):
        global NOT_EQUAL_LENGTH
        NOT_EQUAL_LENGTH.append('{}/{}'.format(species_folder, file_name))
    elif not group <= len(lengths) <= species:
        global BROKEN_SPECIES
        # BROKEN_SPECIES.append('{}/{}'.format(species_folder, file_number))
        if not BROKEN_SPECIES.get(species_folder):
            BROKEN_SPECIES[species_folder] = list()
        BROKEN_SPECIES[species_folder].append(file_name)
    elif lengths[0] % 3 != 0:
        global NOT_MULTIPLE_OF_THREE
        NOT_MULTIPLE_OF_THREE.append('{}/{}'.format(species_folder, file_name))
    else:
        logging.warning("seq lengths are not equal or wrong number of species: {}".format(str(len(lengths))))
        # BROKEN_FILES.append('{}/{}'.format(species_folder, file_number))
        if not BROKEN_LENGTH_FILES.get(species_folder):
            BROKEN_LENGTH_FILES[species_folder] = list()
        BROKEN_LENGTH_FILES[species_folder].append(file_name)


def write_target_phy_file(line_edited, target_file):
    for chunk in chunks(line_edited, 60):
        target_file.write(chunk)


def phylip2paml(folder_out, species_folder, source_file_name, species, group):
    """ converting philip-sequential to specific philip format for paml """
    # file_number = re.search(r'(\d+)\.', source_file_name).group(1)
    file_gene_name = source_file_name.split('.')[0]
    target_file_path = os.path.join(folder_out, species_folder, file_gene_name, '{}.{}'.format(file_gene_name,
                                                                                               "phy"))  # sort for
    # additional folders for paml
    source_file_path = os.path.join(folder_out, species_folder, source_file_name)
    if not os.path.isdir(os.path.join(folder_out, species_folder, file_gene_name)):
        os.makedirs(os.path.join(folder_out, species_folder, file_gene_name))
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

    check_lengths(lengths, species_folder, file_gene_name, species, group)
    logging.info('changing for paml and SWAMP .phy format file {} has been recorded'.format(target_file_path))


def replace_broken_files_and_write_table(directory_out):
    broken_length_folder = os.path.join(directory_out, "broken_length_files")
    broken_species_folder = os.path.join(directory_out, "broken_species_files")
    if BROKEN_LENGTH_FILES:
        for folder, files in BROKEN_LENGTH_FILES.items():
            for f in files:
                shutil.move(os.path.join(directory_out, folder, f), os.path.join(broken_length_folder, folder, f))
        logging.warning("BROKEN_LENGTH_FILES files moves to to {}".format(broken_length_folder))
    if BROKEN_SPECIES:
        for folder, files in BROKEN_SPECIES.items():
            for f in files:
                shutil.move(os.path.join(directory_out, folder, f), os.path.join(broken_species_folder, folder, f))
            # os.remove(os.path.join(directory_out, folder)) # TODO: shutil does not delete source, just leave empty
        logging.warning("BROKEN_SPECIES files moves to to {}".format(broken_species_folder))


def main(folder_in, folder_order, folder_out, species, group):
    for input_tuple in get_infile_and_order(folder_in, folder_order):
        if input_tuple:
            species_folder, infile, order_string = input_tuple
            fasta2phylip(species_folder, infile, order_string, folder_in, folder_out)
    for species_folder, phylip_file in parse_phylip_dir(folder_out):
        phylip2paml(folder_out, species_folder, phylip_file, species, group)
        seq_philip_file = os.path.join(folder_out, species_folder, phylip_file)
        os.remove(seq_philip_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the input folder with fasta files sorted by separated folders', nargs='?',
                        required=True)
    parser.add_argument('--order', help='Path to the folder with .order files for each folder in input', nargs='?',
                        required=True)
    parser.add_argument('--o', help='Path to the output folder with result philip files, if it does not exist,'
                                    ' it will be created automatically', nargs='?', required=True)
    parser.add_argument('--species', help='Number of species', nargs='?', required=True)
    parser.add_argument('--group', help='Minimal size of species group', nargs='?', required=True)
    args = parser.parse_args()
    out_dir = args.o
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    try:
        main(args.i, args.order, out_dir, int(args.species), int(args.group))
        logging.warning("BROKEN_FILES {}".format(len(BROKEN_LENGTH_FILES)))
        if BROKEN_LENGTH_FILES or BROKEN_SPECIES:
            replace_broken_files_and_write_table(out_dir)
        logging.warning("NOT_EQUAL_LENGTH {}".format(len(NOT_EQUAL_LENGTH)))
        logging.warning("BROKEN_SPECIES groups {}".format(len(BROKEN_SPECIES)))
        logging.warning("NOT_MULTIPLE_OF_THREE {}".format(len(NOT_MULTIPLE_OF_THREE)))
        resulting_file = os.path.join(out_dir, 'fasta2paml_ordering_summary.xlsx')
        writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
        df = pd.DataFrame({'Gene name': NOT_EQUAL_LENGTH})
        df.to_excel(writer, sheet_name='NOT_EQUAL_LENGTH', index=False)
        df = pd.DataFrame({'Gene name': NOT_MULTIPLE_OF_THREE})
        df.to_excel(writer, sheet_name='NOT_MULTIPLE_OF_THREE', index=False)
        df = pd.DataFrame({'Gene name': BROKEN_LENGTH_FILES})
        df.to_excel(writer, sheet_name='BROKEN_FILES', index=False)
        df = pd.DataFrame.from_dict(BROKEN_SPECIES)
        df.to_excel(writer, sheet_name='BROKEN_SPECIES', index=False)
        writer.save()
        logging.info("Summary has been written to {}".format(resulting_file))
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
    logging.info("The work has been completed")
