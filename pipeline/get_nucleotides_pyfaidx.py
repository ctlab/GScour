#!/usr/bin/env python
import os
from pyfaidx import Fasta
import argparse
import pandas as pd


def parse_dir(in_dir):
    for f in os.scandir(in_dir):
        if f.is_file() and f.name.split('.')[-1] == 'fna':
            yield f


def main(target_genes_file, input_fasta_dir, output_dir):
    df = pd.io.excel.read_excel(target_genes_file)
    target_genes = list(df['gene name'])
    print("target genes of len={}".format(len(target_genes)))
    for f in parse_dir(input_fasta_dir):
        count = 0
        genes_longest = Fasta(os.path.join(input_fasta_dir, f.name), duplicate_action="longest")
        for gene in target_genes:
            if gene in genes_longest.keys() or genes_longest[gene]:
                with open(os.path.join(output_dir, "{}".format(f.name)), 'a') as fas:
                    fas.write('>' + gene)
                    fas.write(genes_longest[gene])  # or str(seqFile)
                    count += 1
        print("count for {}={}".format(f.name, count))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--t', help='Full path to the .xlsx with target names in column \'gene name\'',
                        nargs='?', required=True)
    parser.add_argument('--i', help='Full path to the fasta files folder',
                        nargs='?', required=True)
    parser.add_argument('--o', help='Full path to the folder with output fasta files', nargs='?',
                        required=True)
    args = parser.parse_args()
    if not os.path.isdir(args.o):
        os.makedirs(args.o)

    main(args.t, args.i, args.o)

