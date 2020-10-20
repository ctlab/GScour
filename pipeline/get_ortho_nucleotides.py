#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

FILE_NUMBER_OF_SEQS = dict()
BROKEN_SPECIES = list()
BROKEN_LENGTH = list()
LOG_FILE = "get_ortho_nuc.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def get_and_write_nucleotide_seq(gb_file, protein_ids, csv_columns, directory_out, species_numerating):
    anti_repeat_store = dict()
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("protein_id") in protein_ids:
                    protein_id = feature.qualifiers.get("protein_id")[0]
                    nucleotide_seq = feature.location.extract(record.seq)
                    seq_length = len(nucleotide_seq)
                    index_count = np.where(protein_ids == protein_id)[0][0]
                    file_out_number = str(index_count + 1)
                    """check if length of nucleotide sequence equal 3*number of proteins"""
                    if protein_id in csv_columns['Protein product'].values:
                        protein_length = csv_columns.loc[csv_columns['Protein product'] == protein_id, 'Length'].values[0]
                        if not seq_length == 3*protein_length:
                            logging.warning("length {} of nucleotide_seq {} not equal proteins length"
                                            "with protein_id {}".
                                            format(seq.id, seq_length, protein_id))
                            BROKEN_LENGTH.append(file_out_number)
                            continue

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

                    gene = feature.qualifiers.get('gene')[0]
                    logging.info("have found gene {} corresponding protein_id {} with nucleotide seq\n{}\n of "
                                 "length".format(gene, protein_id, nucleotide_seq, seq_length))
                    seq = SeqRecord(nucleotide_seq, id=species_numerating, description="")

                    with open(os.path.join(directory_out, file_out_number + ".fna"), "a") as ofile:
                        SeqIO.write(seq, ofile, "fasta")
                        logging.info("wrote nucleotide_seq of length {} in file number {}".
                                     format(seq_length, file_out_number))
                        if not FILE_NUMBER_OF_SEQS.get(file_out_number):
                            FILE_NUMBER_OF_SEQS[file_out_number] = 0
                        FILE_NUMBER_OF_SEQS[file_out_number] += 1

                    with open(os.path.join(directory_out, file_out_number + ".log"), "a") as ofile:
                        ofile.write("{} - {} - {} - {} - {}\n".format(gene, protein_id, seq_length,
                                                                      file_out_number, species_numerating))
                        logging.info("wrote seq of gene {} corresponding {} in file number {}".format(gene, protein_id,
                                                                                                      file_out_number))


def main(orthodata_filepath, annotation_gbff, annotation_csv, species, directory_out):
    if not os.path.isdir(directory_out):
        os.makedirs(directory_out)
    ortho_data = pd.read_csv(orthodata_filepath, sep='\t', usecols=range(3, 3 + species))
    for column_number, _ in enumerate(ortho_data.columns):
        species_numerating = str(column_number + 1)
        annotation_gbff_path = os.path.join(annotation_gbff, "{}.{}".format(species_numerating, 'gbff'))
        annotation_csv_path = os.path.join(annotation_csv, "{}.{}".format(str(species_numerating), 'csv'))
        logging.info('Working with annotation files {}\n{}'.format(annotation_gbff_path, annotation_csv_path))
        logging.info('Working with species number {}'.format(str(species_numerating)))
        df = pd.read_csv(annotation_gbff_path, sep=',', names=['Protein product', 'Length'], usecols=[8, 9])
        get_and_write_nucleotide_seq(annotation_gbff_path, ortho_data.iloc[:, column_number].values, df, directory_out,
                                     species_numerating)
    for file, species in FILE_NUMBER_OF_SEQS.items():
        if species != 5:
            BROKEN_SPECIES.append(file)


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
        if BROKEN_SPECIES:
            logging.warning("BROKEN_SPECIES: {}".format(BROKEN_SPECIES))
        if BROKEN_LENGTH:
            logging.warning("BROKEN_LENGTH: {}".format(BROKEN_LENGTH))
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
