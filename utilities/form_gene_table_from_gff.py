#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from BCBio import GFF
import pandas as pd
import os

""" script to make list of genes from .gff file """

gff_genes = dict()


def parse_dir(input_dir):
    if os.path.isdir(input_dir):
        for infile in os.scandir(input_dir):
            if infile.name.split('.')[-1] == 'gff':
                yield os.path.join(input_dir, infile.name), infile.name.split('.')[0]


def main(gff_folder, out_file):
    for infile_gff_path, file_name in parse_dir(gff_folder):
        with open(infile_gff_path, 'r') as f:
            for rec in GFF.parse(f):
                for feature in rec.features:
                    if feature.type == 'gene':
                        if feature.qualifiers.get('gene'):
                            gene = feature.qualifiers.get('gene')[0]
                        if not gff_genes.get(file_name):
                            gff_genes[file_name] = set()
                        gff_genes[file_name].add(gene)

    summary_df = pd.DataFrame()
    for k, v in gff_genes.items():
        df = pd.DataFrame({k: sorted(list(v))})
        summary_df = pd.concat([summary_df, df], axis=1, sort=True)
    writer = pd.ExcelWriter(out_file, engine='openpyxl')
    summary_df.to_excel(writer, sheet_name='genes from .gff')
    writer.save()
    print("done, table has been recorded to {}".format(out_file))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder with .gff files', nargs='?', required=True)
    args = parser.parse_args()
    gff_folder_path = args.i
    out_file_xlsx = os.path.join(gff_folder_path, "gene_table_from_gff.xlsx")
    print("--i option - path to the folder with .gff files is {}\n"
          "path to the output table with genes {}".format(gff_folder_path, out_file_xlsx))
    try:
        main(gff_folder_path, out_file_xlsx)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
