#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pandas import Index

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

DEFAULT_DATA_PATH = os.getcwd()
DEFAULT_FOLDER_OUT = os.path.join(os.getcwd(), 'ortho_fna/')


def write_fasta_files(initfna_filepath, ortho_data, column_number, gene_id, start, stop, protein_id, directory_out):
    for record in SeqIO.parse(initfna_filepath, "fasta"):
        if record.id == gene_id:
            logging.info("gene_id {0} has been found in file {1}".format(record.id, initfna_filepath))
            index_count = Index(ortho_data.iloc[:, column_number]).get_loc(protein_id)
            seq = SeqRecord(record.seq[start:stop + 1], id=record.id, description=record.description,
                            annotations=record.annotations)
            with open(os.path.join(directory_out, str(index_count) + ".fna"), "a") as ofile:
                logging.info("writing file {}".format(ofile.name))
                SeqIO.write(seq, ofile, "fasta")
            logging.info("wrote gene_id {} in file number {}".format(gene_id, index_count))
            break


def main(initfna_folder, annotation_folder, orthodata_filepath, species, directory_out):
    if not os.path.isdir(directory_out):
        os.makedirs(directory_out)
    ortho_data = pd.read_csv(orthodata_filepath, sep='\t', usecols=range(3, 3 + int(species)))
    for column_number, file_faa in enumerate(ortho_data.columns.values):
        species_numerating = column_number + 1
        initfna_filepath = os.path.join(initfna_folder, "{}.{}".format(str(species_numerating),'fna'))
        annotation_filepath = os.path.join(annotation_folder, "{}.{}".format(str(species_numerating),'csv'))
        logging.info('Working with species number {}'.format(str(species_numerating)))
        df = pd.read_csv(annotation_filepath, sep=',', usecols=[1, 2, 3, 8])
        gene_counter = 0
        anti_repeat_store = dict()
        for protein_id in ortho_data.iloc[:, column_number].values:
            if protein_id in df['Protein product'].values:
                if not anti_repeat_store.get(protein_id):
                    anti_repeat_store[protein_id] = 1
                else:
                    continue
                gene_id = df.loc[df['Protein product'] == protein_id, 'Accession'].values[0]
                start = df.loc[df['Protein product'] == protein_id, 'Start'].values[0]
                stop = df.loc[df['Protein product'] == protein_id, 'Stop'].values[0]
                logging.info("find_gene_protein_ids {}, {}".format(gene_id, protein_id))
                write_fasta_files(initfna_filepath, ortho_data, column_number, gene_id,
                                  start, stop, protein_id, directory_out)
                gene_counter += 1
        logging.info("{} of nucleotide sequences were recorded from init fna file".
                     format(gene_counter, initfna_filepath))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--initfna', help='Path to the folder with genome init files .fna from '
                                          'www.ncbi.nlm.nih.gov/genome/', nargs='?',
                        default=DEFAULT_DATA_PATH)
    parser.add_argument('--annotation', help='Path to the folder with annotation .csv files from '
                                             'www.ncbi.nlm.nih.gov/genome/', nargs='?',
                        default=DEFAULT_DATA_PATH)
    parser.add_argument('--ortho', help='Path to the single_copy_orthologs.tsv', nargs='?',
                        default=DEFAULT_DATA_PATH)
    parser.add_argument('--species', help='Number of species', nargs='?')
    parser.add_argument('--out', help='Path to the folder for result write out', nargs='?',
                        default=DEFAULT_FOLDER_OUT)
    args = parser.parse_args()

    try:
        main(args.initfna, args.annotation, args.ortho, args.species, args.out)
    except:
        logging.exception("Unexpected error")
