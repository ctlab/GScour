#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

WRITTEN_FILES = dict()
BROKEN_SPECIES = list()
BROKEN_LENGTH = list()
NUMBER_OF_NEED_TO_BE_WRITTEN = 0
LOG_FILE = "get_ortho_nuc.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def get_and_write_nucleotide_seq(gb_file, ortho_protein_ids, csv_columns, directory_out, species_numerating):
    global WRITTEN_FILES
    global BROKEN_LENGTH
    anti_repeat_store = dict()
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("protein_id") in ortho_protein_ids:
                    protein_id = feature.qualifiers.get("protein_id")[0]
                    nucleotide_seq = feature.location.extract(record.seq)
                    seq_length = len(nucleotide_seq)
                    index_count = np.where(ortho_protein_ids == protein_id)[0][0]
                    file_out_number = str(index_count + 1)
                    seq = SeqRecord(nucleotide_seq, id=species_numerating, description="")
                    """check if length of nucleotide sequence equal 3*number of proteins"""
                    if protein_id in csv_columns['Protein product'].values:
                        protein_length = csv_columns.loc[csv_columns['Protein product'] == protein_id, 'Length'].values[
                            0]
                        if not ((seq_length - 3) == 3 * protein_length):  # except stop-codon
                            logging.warning("length {} of nucleotide_seq {} not equal proteins length {}*3"
                                            "with protein_id {}, species {}".
                                            format(str(seq_length - 3), seq.id, protein_length, protein_id,
                                                   species_numerating))
                            if file_out_number not in BROKEN_LENGTH:
                                BROKEN_LENGTH.append(file_out_number)
                            continue
                    else:
                        logging.warning("No inspections for protein length were carried out")

                    if not len(seq) % 3 == 0:
                        if file_out_number not in BROKEN_LENGTH:
                            BROKEN_LENGTH.append(file_out_number)
                        logging.warning("length {} not multiple of three in file number {}".
                                        format(len(seq), file_out_number))
                        continue
                    else:
                        logging.info("nucleotide_seq of length {} multiple of three for file number {}".
                                     format(str(seq_length - 3), file_out_number))
                    logging.info("LENGTH CHECK: seq {}, seq in object Seq {}".format(seq_length, len(seq)))

                    """ check duplicates"""
                    if not anti_repeat_store.get(protein_id):
                        anti_repeat_store[protein_id] = nucleotide_seq
                    else:
                        logging.warning("REPEAT file {} protein_id {}".format(file_out_number, protein_id))
                        if anti_repeat_store[protein_id] == nucleotide_seq:
                            logging.warning("REPEAT file {} protein_id {} nucleotide_seqs is equal".format(
                                file_out_number, protein_id))
                        else:
                            logging.warning("REPEAT file {} protein_id {} nucleotide_seqs is NOT equal".format(
                                file_out_number, protein_id))
                        continue

                    if file_out_number in BROKEN_LENGTH:
                        continue

                    gene = feature.qualifiers.get('gene')[0]
                    logging.info("have found gene {} corresponding protein_id {} with nucleotide seq\n{}\n of "
                                 "length".format(gene, protein_id, nucleotide_seq, seq_length))

                    with open(os.path.join(directory_out, file_out_number + ".fna"), "a") as ofile:
                        SeqIO.write(seq, ofile, "fasta")
                        if not WRITTEN_FILES.get(file_out_number):
                            WRITTEN_FILES[file_out_number] = 0
                        WRITTEN_FILES[file_out_number] += 1

                    log_file = os.path.join(directory_out, file_out_number + ".log")
                    with open(log_file, "a") as ofile:
                        filesize = os.path.getsize(log_file)
                        if not filesize == 0:
                            ofile.write("gene_name - protein_id - seq_length - file number - species_numerating\n")
                        ofile.write("{} - {} - {} - {} - {}\n".format(gene, protein_id, seq_length,
                                                                      file_out_number, species_numerating))
                        logging.info("wrote seq of gene {} corresponding {} in file number {}".format(gene, protein_id,
                                                                                                      file_out_number))


def conv_int(val):
    if not val:
        return 0
    if val == 'Length':
        return 'Length'
    try:
        return np.int(val)
    except:
        return np.int(0)


def conv_string(val):
    if not val:
        return ""
    if val == 'Protein product':
        return 'Protein product'
    try:
        return np.str(val)
    except:
        return np.str("")


def replace_broken_files(directory_out):
    broken_folder = "broken_species_files"
    os.makedirs(broken_folder)
    for file_number in BROKEN_SPECIES:
        os.replace(os.path.join(directory_out, file_number + ".fna"), os.path.join(broken_folder, file_number + ".fna"))
        os.replace(os.path.join(directory_out, file_number + ".log"), os.path.join(broken_folder, file_number + ".log"))


def main(orthodata_filepath, annotation_gbff, annotation_csv, species, directory_out):
    global BROKEN_SPECIES
    global NUMBER_OF_NEED_TO_BE_WRITTEN
    if not os.path.isdir(directory_out):
        os.makedirs(directory_out)
    ortho_data = pd.read_csv(orthodata_filepath, sep='\t', usecols=range(3, 3 + species))
    for column_number, _ in enumerate(ortho_data.columns):
        species_numerating = str(column_number + 1)
        annotation_gbff_path = os.path.join(annotation_gbff, "{}.{}".format(species_numerating, 'gbff'))
        annotation_csv_path = os.path.join(annotation_csv, "{}.{}".format(species_numerating, 'csv'))
        logging.info('Working with annotation files {}\n{}'.format(annotation_gbff_path, annotation_csv_path))
        logging.info('Working with species number {}'.format(species_numerating))
        df = pd.read_csv(annotation_csv_path, sep=',', names=['Protein product', 'Length'], usecols=[8, 9],
                         converters={'Protein product': conv_string, 'Length': conv_int}, low_memory=False)
        ortho_protein_ids = ortho_data.iloc[:, column_number].values

        NUMBER_OF_NEED_TO_BE_WRITTEN = len(ortho_protein_ids)
        logging.info("NUMBER_OF_NEED_TO_BE_WRITTEN: {}".format(NUMBER_OF_NEED_TO_BE_WRITTEN))

        get_and_write_nucleotide_seq(annotation_gbff_path,ortho_protein_ids , df, directory_out,
                                     species_numerating)

        if WRITTEN_FILES:
            logging.info("WRITTEN_FILES {}:\n{}".format(len(WRITTEN_FILES), repr(WRITTEN_FILES)))
        if BROKEN_LENGTH:
            logging.warning("BROKEN_LENGTH: {}".format(BROKEN_LENGTH))
    for file, written_species in WRITTEN_FILES.items():
        if written_species != species:
            if file not in BROKEN_SPECIES:
                BROKEN_SPECIES.append(file)
    if BROKEN_SPECIES:
        logging.warning("BROKEN_SPECIES: {}".format(BROKEN_SPECIES))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Path to the single_copy_orthologs.tsv', nargs='?')
    parser.add_argument('--gbff', help='Path to the folder with annotation .gbff files from '
                                       'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--csv', help='Path to the folder with annotation .csv files from '
                                      'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--species', help='Number of species', nargs='?')
    parser.add_argument('--out', help='Path to the folder for result write out', nargs='?')
    args = parser.parse_args()

    try:
        main(args.ortho, args.gbff, args.csv, int(args.species), args.out)
        written_files_number = len(WRITTEN_FILES)
        delta = NUMBER_OF_NEED_TO_BE_WRITTEN - written_files_number
        if delta == 0:
            logging.info("All files are written")
        else:
            logging.info("NUMBER_OF_NEED_TO_BE_WRITTEN = {},  WRITTEN_FILES = {}, where {} in BROKEN_SPECIES list: {}"
                         .format(NUMBER_OF_NEED_TO_BE_WRITTEN, written_files_number, len(BROKEN_SPECIES),
                                 repr(BROKEN_SPECIES)))
            if BROKEN_SPECIES:
                replace_broken_files(args.out)
                residue = written_files_number - len(BROKEN_SPECIES)
                logging.info("removed broken species files into folder 'broken_species_files' in cwd,"
                             "please check out folder for .fna files number: {}".format(residue))

    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
