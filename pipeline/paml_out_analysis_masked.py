#!/usr/bin/env python
import argparse
import sys
import logging
from scipy import stats
import os
import logging
import re
import traceback

LN_NP_PATTERN = re.compile(r"lnL\(ntime:\s+\d+\s+np:\s+(\d+)\):\s+(-\d+\.\d+)")
# POSITIVE_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
POS_SITES_STRING = re.compile(r"Positive\ssites\sfor\sforeground\slineages\sProb\(w>1\):")
POSITION_ACID_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")


# LOG_FILE = "analyse_paml_out_masked.log"


def get_ln_np(infile):
    with open(infile, "r") as f:
        pos_sites = list()
        np = 0
        ln = 0.
        for line in f:
            if re.search(LN_NP_PATTERN, line):
                np = int(re.search(LN_NP_PATTERN, line).group(1))
                ln = float(re.search(LN_NP_PATTERN, line).group(2))
            if re.search(POS_SITES_STRING, line):
                pattern = re.search(POSITION_ACID_PATTERN, f.readline())
                while pattern:
                    position = int(pattern.group(1))
                    acid = pattern.group(2)
                    probability = float(pattern.group(3))
                    pos_sites.append([position, acid, probability])
                    pattern = re.search(POSITION_ACID_PATTERN, f.readline())
    return np, ln, pos_sites


def calc_p_value(np0, ln0, np1, ln1):
    ln_val = 2 * (float(ln1) - (float(ln0)))
    n = np1 - np0
    p_val = 1 - stats.chi2.cdf(ln_val, n)
    return p_val


def count_sites(in_folder, species_folder, item, child_logger):
    broken_paml_outs = list()
    no_significance = 0
    positive_sites_number = 0
    positive_sites_dict = dict()
    positive_genes = list()

    np0, ln0, np1, ln1 = 0, 0., 0, 0.
    pos_sites = []
    if os.path.isdir(item):
        item_folder_name = item.name
        for infile in os.listdir(item):
            if infile.endswith("_null1_masked.out"):
                np0, ln0, _ = get_ln_np(os.path.join(in_folder, species_folder.name,
                                                     item_folder_name, infile))
            if infile.endswith("_alter1_masked.out"):
                np1, ln1, pos_sites = get_ln_np(os.path.join(in_folder, species_folder.name,
                                                             item_folder_name, infile))
        if all([np0, np1, ln0, ln1]):
            p_val = calc_p_value(np0, ln0, np1, ln1)
            if p_val and p_val < 0.05:
                number_pos = len(pos_sites)
                child_logger.info("{} number of positive sites {}".format(item_folder_name, number_pos))
                positive_sites_number += number_pos
            for sites in pos_sites:
                pos, acid, probability = [sites[i] for i in range(3)]
                child_logger.info("{} positive sites : position, acid, probability : {}, {}, {}".
                                  format(item_folder_name, pos, acid, probability))
                if item_folder_name not in positive_genes:
                    positive_genes.append(item_folder_name)
                if not positive_sites_dict.get(item_folder_name):
                    positive_sites_dict[item_folder_name] = pos_sites
            else:
                child_logger.info("{} no significance, p-value {}".format(item_folder_name, p_val))
                no_significance += 1
        else:
            child_logger.warning("lack of params {}: np0 {}, ln0 {}, np1 {}, ln1 {}, pos_sites {}".format(
                item_folder_name, np0, ln0, np1, ln1, pos_sites))
            broken_paml_outs.append(item_folder_name)
    return broken_paml_outs, no_significance, positive_sites_number, positive_sites_dict, positive_genes


def main(in_folder, ortho_logs):
    for species_folder in os.scandir(in_folder):
        broken_paml_outs = list()
        no_significance = 0
        positive_sites_number = 0
        positive_sites_dict = dict()
        positive_genes = list()
        if os.path.isdir(species_folder):
            result_file = os.path.join(in_folder, species_folder.name, "{}.{}".format(species_folder.name,
                                                                                      "result"))
            child_logger = logging.getLogger(__name__)
            child_logger.addHandler(logging.FileHandler(result_file))
            child_logger.setLevel(10)
            for item in os.scandir(species_folder):
                broken_paml_outs_item, no_significance_item, positive_sites_number_item, positive_sites_dict_item, \
                    positive_genes_item = count_sites(in_folder, species_folder, item, child_logger)
                broken_paml_outs += broken_paml_outs_item
                no_significance += no_significance_item
                positive_sites_number += positive_sites_number_item
                print("positive_sites_dict", positive_sites_dict)
                positive_sites_dict.update(positive_sites_dict_item)
                print("positive_sites_dict_item", positive_sites_dict_item)
                positive_genes += positive_genes_item

            child_logger.warning("broken_paml_outs : {} : {}".format(len(broken_paml_outs), broken_paml_outs))
            child_logger.info("Number of no significance files {}".format(no_significance))
            child_logger.info("Number of positive sites {}".format(positive_sites_number))
            child_logger.info("Number of positive genes {} : {}".format(len(positive_genes), positive_genes))
            child_logger.info(
                    "Positive sites : file : position, acid, probability\n{}".format(repr(positive_sites_dict)))

            if ortho_logs:
                gene_names_dict = get_names_of_genes_under_positive(positive_genes, ortho_logs, child_logger)
                if gene_names_dict:
                    child_logger.info("Gene_name_dict of length {}:\nGenes under positive selection: file: gene, "
                                      "protein:\n{} "
                                      .format(len(gene_names_dict), repr(gene_names_dict)))


def get_names_of_genes_under_positive(genes_under_positive, seq_log_folder, child_logger):
    dict_of_gene_names = dict()
    pattern = ""
    for infile in os.scandir(seq_log_folder):
        file_number = infile.name.split('.')[0]
        if file_number in genes_under_positive and infile.name.endswith("log"):
            with open(os.path.join(seq_log_folder, infile.name), "r") as f:
                for line in f:
                    if re.search(r"-\s\d+$", line):
                        pattern = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)")
                        gene_name = (re.search(pattern, line)).group(1)
                        protein_name = (re.search(pattern, line)).group(2)
                        if not dict_of_gene_names.get(infile):
                            dict_of_gene_names[file_number] = [gene_name, protein_name]
    if not pattern:
        child_logger.info("log pattern not found to get gene_name")
    return dict_of_gene_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='The full path to the folder contains folders with input files for paml',
                        nargs='?')
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?', default='')
    args = parser.parse_args()
    infolder = args.infolder
    log_folder = args.log
    print("Passed args: infolder {}, log folder {}".format(infolder, log_folder))
    try:
        main(infolder, log_folder)
    except BaseException as e:
        print("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))
    print("The work has been completed")
