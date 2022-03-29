#!/usr/bin/env python
import argparse
import os
import re
import pandas as pd
"""
the script takes result table with "summary" sheet, extract values from column "Gene name"
and return table with pairs Gene name-number of file with the appropriate sequence
"""

LOG_PATTERN = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)\s-\s(\d+)\s-\s(\d+)\s-\s(\d)")


def get_gene_name_file_number(seq_log_folder):
    """ get gene name and protein_id which corresponds to the required_species """
    pattern, gene_name, protein_id = "", "", ""
    for infile in os.scandir(seq_log_folder):
        file_number = infile.name.split('.')[0]
        if infile.name.endswith("log"):
            print("open log file {}".format(os.path.join(seq_log_folder, infile.name)))
            with open(os.path.join(seq_log_folder, infile.name), "r") as f:
                i = 0
                for line in f:
                    if i == 0:
                        i += 1
                        continue
                    try:
                        gene_name = (re.search(LOG_PATTERN, line)).group(1)
                        print(gene_name)
                        file_number = (re.search(LOG_PATTERN, line)).group(4)
                        yield gene_name, file_number
                        break
                    except AttributeError:
                        if not gene_name:
                            print("ERROR: Can't parse gene name from log")
                            gene_name = "gene-{}".format(file_number)


def main(input_table_path, out_folder, log_folder):
    gene_names = list()
    file_numbers = list()
    df = pd.io.excel.read_excel(input_table_path, sheet_name='summary')
    coincidence = 0
    for gene_name, file_number in get_gene_name_file_number(log_folder):
        if gene_name in df['Gene name'].tolist():
            coincidence += 1
            gene_names.append(gene_name)
            file_numbers.append(file_number)
    print("coincidence=", coincidence)
    result_sheet = {'Gene name': gene_names, 'File number': file_numbers}
    df = pd.DataFrame(result_sheet, columns=['Gene name', 'File number'])
    result_sheet_path = os.path.join(out_folder, 'gene name_file number.xlsx')
    writer = pd.ExcelWriter(result_sheet_path, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='summary')
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the .xlsx table with summary sheet',
                        nargs='?', required=True)
    parser.add_argument('--o', help='The full path to the output folder for writing results',
                        nargs='?', required=True)
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?', required=True)
    args = parser.parse_args()
    if not os.path.isdir(args.o):
        os.makedirs(args.o)
    main(args.i, args.o, args.log)
    print("done")
