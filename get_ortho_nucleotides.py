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
directory_init_fna = '/scratch/afedorova/acfhp/acfhp_GCF_fna'
directory_out = '/scratch/afedorova/acfhp/acfhp_ortho_fna'
csv_dir = "/scratch/afedorova/acfhp/csv/"
directory_tsv = "/home/alina_grf/progprojects/data_work/"
ortho_data = pd.read_csv(os.path.join(directory_tsv, 'single_copy_orthologs.tsv'), sep='\t', usecols=[3, 4, 5, 6, 7])


def get_names(file_faa):
    if file_faa.split('.')[1] == '39_GRCh38':
        file_fna = file_faa.replace('p13_protein.faa', 'p13_genomic.fna')
        fna_filepath = os.path.join(directory_init_fna, file_fna)
        csv_file = os.path.join(csv_dir, file_fna.replace('fna', 'csv'))
    if file_faa.split('.')[1] == '3_CanFam3':
        file_fna = file_faa.replace('1_protein.faa', '1_genomic.fna')
        fna_filepath = os.path.join(directory_init_fna, file_fna)
        csv_file = os.path.join(csv_dir, file_fna.replace('fna', 'csv'))
    if file_faa.split('.')[1] == '3_Felis_catus_9':
        file_fna = file_faa.replace('0_protein.faa', '0_genomic.fna')
        fna_filepath = os.path.join(directory_init_fna, file_fna)
        csv_file = os.path.join(csv_dir, file_fna.replace('fna', 'csv'))
    if file_faa.split('.')[1] == '1_PanTig1':
        file_fna = file_faa.replace('0_protein.faa', '0_genomic.fna')
        fna_filepath = os.path.join(directory_init_fna, file_fna)
    if file_faa.split('.')[1] == '1_Aci_jub_2_genomic':
        file_fna = file_faa.replace('faa', 'fna')
        fna_filepath = os.path.join(directory_init_fna, file_fna)
        csv_file = os.path.join(csv_dir, file_fna.replace('fna', 'csv'))
    return file_fna, fna_filepath, csv_file


def write_fasta_files(fna_filepath, column_number, gene_id, start, stop, protein_id):
    for record in SeqIO.parse(os.path.join(directory_init_fna, fna_filepath), "fasta"):
        if record.id == gene_id:
            logging.info("gene_id {0} has been found in file {1}".format(record.id, fna_filepath))
            index_count = Index(ortho_data.iloc[:, column_number]).get_loc(protein_id)
            seq = SeqRecord(record.seq[start:stop + 1], id=record.id, description=record.description,
                            annotations=record.annotations)
            with open(os.path.join(directory_out, str(index_count) + ".fna"), "a") as ofile:
                logging.info("writing file {}".format(ofile.name))
                SeqIO.write(seq, ofile, "fasta")
            logging.info("wrote gene_id {} in file number {}".format(gene_id, index_count))
            break


def main():
    for column_number, file_faa in enumerate(ortho_data.columns.values):
        species_numerating = column_number + 1
        file_fna, fna_filepath, csv_file = get_names(file_faa)
        logging.info('Working with file {0}: "{1}"'.format(column_number, csv_file))
        df = pd.read_csv(csv_file, sep=',', usecols=[1, 2, 3, 8])
        gene_counter = 0
        anti_repeat_store = dict()
        """
        if column_number == 0:
            ls = [i for i in range(0, 1)] + [i for i in range(798, 3560)]
            """
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
                write_fasta_files(fna_filepath, column_number, gene_id, start, stop, protein_id)
                gene_counter += 1
        logging.info("gene counter equal {} with file {}".format(gene_counter, csv_file))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--initfna', help='Path to the folder with genome init files .fna from '
                                          'www.ncbi.nlm.nih.gov/genome/', nargs='?',
                        default=DEFAULT_DATA_PATH)
    parser.add_argument('--annotation', help='Path to the annotation .csv file from '
                                             'www.ncbi.nlm.nih.gov/genome/', nargs='?',
                        default=DEFAULT_DATA_PATH)
    parser.add_argument('--ortho', help='Path to the single_copy_orthologs.tsv ', nargs='?',
                        default=DEFAULT_DATA_PATH)
    parser.add_argument('--out', help='Path to the folder for result write out ', nargs='?',
                        default=DEFAULT_FOLDER_OUT)
    args = parser.parse_args()


    if not os.path.isdir(directory_out):
        os.makedirs(directory_out)
    try:
        main()
    except:
        logging.exception("Unexpected error")
