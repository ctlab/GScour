#!/usr/bin/env python
import pandas as pd
import argparse
import re
import os
import logging
"""
this script checks logs for absent items from .xlsx table
with one column - gene names
"""


def get_item_by_gene_name(gene_name, log_folder_path, target_species, logger):
    """ get item by gene name of the required_species """
    for infile in os.scandir(log_folder_path):
        file_number = infile.name.split('.')[0]
        if infile.name.endswith("log"):
            # logger.info("open log file {}".format(os.path.join(log_folder_path, infile.name)))
            with open(os.path.join(log_folder_path, infile.name), "r") as f:
                for line in f:
                    if re.search(r"-\s[{}]+$".format(target_species), line):
                        pattern = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)")
                        gene_name_log = (re.search(pattern, line)).group(1)
                        protein_id_log = (re.search(pattern, line)).group(2)
                        if gene_name == gene_name_log:
                            logger.info('found {} in {}, protein_id {}'.format(gene_name, file_number, protein_id_log))
                            return file_number, protein_id_log
    return gene_name


def main(excel_file_path, log_folder_path, target_species, logger):
    absent_genes = []
    genes_items = []
    df = pd.io.excel.read_excel(excel_file_path)
    child_logger.info("genes to check of length {}".format(len(df)))
    for gene_name_tuple in df.itertuples():
        gene_name = gene_name_tuple[1]
        returned = get_item_by_gene_name(gene_name, log_folder_path, target_species, logger)
        if returned == gene_name:
            absent_genes.append(gene_name)
        else:
            genes_items.append((gene_name, returned[1], returned[0]))
    logger.info("absent genes in ortho log files {}:\n{}".format(len(absent_genes), absent_genes))
    logger.info("genes which was found in items {}:\n{}".format(len(genes_items), genes_items))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the .xlsx file with one column of gene names to check', nargs='?')
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?')
    parser.add_argument('--required', help='Number of required (single target) species for analysis', nargs='?')
    args = parser.parse_args()
    in_file = args.i
    log_folder = args.log
    required_species = args.required
    log_file = os.path.join(os.path.split(in_file)[0], "{}.{}".format("check_logs", "log"))
    child_logger = logging.getLogger('__main__.' + in_file)
    child_logger.addHandler(logging.FileHandler(log_file))
    child_logger.setLevel(10)
    child_logger.info("Passed args: input directory {}, log folder {}, required species {}".format(in_file, log_folder,
                                                                                                   required_species))
    child_logger.info("Log file: {}".format(log_file))

    try:
        main(in_file, log_folder, required_species, child_logger)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
    print("The work has been completed")

