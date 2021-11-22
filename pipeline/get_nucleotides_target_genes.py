#!/usr/bin/env python
import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_dir(in_dir):
    for f in os.scandir(in_dir):
        if f.is_file() and f.name.split('.')[-1] == 'gbff':
            yield f


def get_max_seq_from_gbff(gb_file, target_genes_series):
    lengths = dict()
    gene_seq = dict()
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get('gene')[0] in list(target_genes_series):
                    gene_name = feature.qualifiers.get("gene")[0]
                    protein_id = feature.qualifiers.get("protein_id")[0]
                    protein_translation_length = len(feature.qualifiers.get("translation")[0])  # int
                    nucleotide_seq = feature.location.extract(record.seq)
                    index_count = str(np.where(target_genes_series == gene_name)[0][0])
                    if not lengths.get(gene_name):
                        lengths[gene_name] = {protein_translation_length: {'protein_id': protein_id,
                                              'seq': nucleotide_seq, 'idx': index_count,
                                                                           'lengths': protein_translation_length}}
                    else:
                        lengths[gene_name].update({
                            protein_translation_length: {
                                'protein_id': protein_id,
                                'seq': nucleotide_seq, 'idx': index_count, 'lengths': protein_translation_length
                                }
                            })
    for k, v_dict in lengths.items():
        max_key = max(v_dict.keys())
        gene_seq[k] = lengths[k][max_key]
    return gene_seq


def write_fasta_and_log_files(directory_out, gene_seq_dict, species_numerating):
    count = 0
    for k, v in gene_seq_dict.items():
        seq_record = SeqRecord(v['seq'], id=species_numerating, description="")
        with open(os.path.join(directory_out, '{}.{}'.format(v['idx'], 'fna')), "a") as o_file:
            SeqIO.write(seq_record, o_file, "fasta")
            count += 1
        log_file = os.path.join(directory_out, '{}.{}'.format(v['idx'], 'log'))
        with open(log_file, "a") as o_file:
            file_size = os.path.getsize(log_file)
            if file_size == 0:
                o_file.write("gene_name - protein_id - seq_length - file number - species_numerating\n")
            o_file.write("{} - {} - {} - {} - {}\n".format(k, v['protein_id'], v['lengths'],
                                                          v['idx'], species_numerating))
            print("Wrote gene {}, protein_id {}, seq_length {}, file_out_number {}, species_numerating {}"
                  .format(k, v['protein_id'], v['lengths'], v['idx'], species_numerating))
    print("Sum of genes wrote = {}".format(count))


def main(target_genes_file, input_gbff_dir, output_dir):
    df = pd.io.excel.read_excel(target_genes_file)
    print("target genes of len={}".format(len(df['gene name'])))
    for gbff_file in parse_dir(input_gbff_dir):
        gene_seq = get_max_seq_from_gbff(gbff_file, df['gene name'])
        write_fasta_and_log_files(output_dir, gene_seq, gbff_file.name.split('.')[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--t', help='Full path to the .xlsx with target names in column \'gene name\'',
                        nargs='?', required=True)
    parser.add_argument('--gbff', help='Path to the folder with annotation .gbff files from '
                                       'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--o', help='Full path to the folder with output fasta files', nargs='?',
                        required=True)
    args = parser.parse_args()
    if not os.path.isdir(args.o):
        os.makedirs(args.o)

    main(args.t, args.gbff, args.o)
