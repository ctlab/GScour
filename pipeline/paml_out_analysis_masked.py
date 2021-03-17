#!/usr/bin/env python
import argparse
import sys
from scipy import stats
import os
import logging
import re
import traceback

LN_NP_PATTERN = re.compile(r"lnL\(ntime:\s+\d+\s+np:\s+(\d+)\):\s+(-\d+\.\d+)")
# POSITIVE_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
POS_SITES_STRING = re.compile(r"Positive\ssites\sfor\sforeground\slineages\sProb\(w>1\):")
POSITION_ACID_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
common_pos_gene_dict = dict()


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
    broken_paml_outs_item = list()
    no_significance_item = 0
    positive_sites_number_item = 0
    items_of_pos_site_dict_item = dict()

    np0, ln0, np1, ln1 = 0, 0., 0, 0.
    pos_sites = []
    if os.path.isdir(item):
        item_folder_name = item.name
        item_id = "{}/{}".format(species_folder.name, item_folder_name)
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
                child_logger.info("{} p-value {}: number of positive sites {}".format(item_id, p_val, number_pos))
                positive_sites_number_item += number_pos
                for sites in pos_sites:
                    pos, acid, probability = [sites[i] for i in range(3)]
                    child_logger.info("{} positive sites : position, acid, probability : {}, {}, {}".
                                      format(item_id, pos, acid, probability))

                if not items_of_pos_site_dict_item.get(item_folder_name):
                    items_of_pos_site_dict_item[item_folder_name] = pos_sites
            else:
                child_logger.info("{} no significance, p-value {}".format(item_id, p_val))
                no_significance_item += 1
        else:
            child_logger.warning("{} lack of params: np0 {}, ln0 {}, np1 {}, ln1 {}, pos_sites {}".format(
                item_id, np0, ln0, np1, ln1, pos_sites))
            broken_paml_outs_item.append(item_folder_name)
    return broken_paml_outs_item, no_significance_item, positive_sites_number_item, items_of_pos_site_dict_item


def main(in_folder, ortho_logs, required_species):
    global common_pos_gene_dict
    for species_folder in os.scandir(in_folder):
        broken_paml_outs = list()
        no_significance = 0
        positive_sites_number = 0
        items_of_pos_site_dict = dict()

        if os.path.isdir(species_folder):
            result_file = os.path.join(in_folder, species_folder.name, "{}.{}".format(species_folder.name,
                                                                                      "result"))
            child_logger = logging.getLogger('__main__.' + species_folder.name)
            child_logger.addHandler(logging.FileHandler(result_file))
            child_logger.setLevel(10)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    broken_paml_outs_item, no_significance_item, positive_sites_number_item, \
                        items_of_pos_site_dict_item = count_sites(in_folder, species_folder, item, child_logger)

                    broken_paml_outs += broken_paml_outs_item
                    no_significance += no_significance_item
                    positive_sites_number += positive_sites_number_item
                    items_of_pos_site_dict.update(items_of_pos_site_dict_item)

            child_logger.warning("Species folder {}: broken_paml_outs : {} : {}".
                                 format(species_folder.name, len(broken_paml_outs), broken_paml_outs))
            child_logger.info("Species folder {}: number of no significance files {}".format(species_folder.name,
                                                                                             no_significance))
            child_logger.info("Species folder {}: number of positive sites {}".format(species_folder.name,
                                                                                      positive_sites_number))
            child_logger.info(
                "Species folder {}: items of positive sites {}: file : position, acid, probability\n{}".format(
                    species_folder.name, len(items_of_pos_site_dict), repr(items_of_pos_site_dict)))
        if os.path.isdir(species_folder):
            if ortho_logs:
                if items_of_pos_site_dict:
                    gene_names_dict = get_names_of_genes_under_positive(items_of_pos_site_dict, ortho_logs,
                                                                        required_species, child_logger)
                    if gene_names_dict:
                        child_logger.info("Gene_name_dict of length {}:\nGenes under positive selection: gene: "
                                          "protein, "
                                          "item:\n{} "
                                          .format(len(gene_names_dict), repr(gene_names_dict)))

                        for key in gene_names_dict.keys():
                            if key not in common_pos_gene_dict.keys():
                                common_pos_gene_dict[key] = list()
                            common_pos_gene_dict[key] += gene_names_dict[key]
                else:
                    child_logger.info("No items of positive sites")
            else:
                child_logger.info("Path to the log folder was not provided")


def get_names_of_genes_under_positive(items_of_pos_site_dict, seq_log_folder, required_species, child_logger):
    """ get gene name and protein_id wich corresponds to the required_species """
    dict_of_gene_names = dict()
    pattern = ""
    for infile in os.scandir(seq_log_folder):
        file_number = infile.name.split('.')[0]
        if file_number in items_of_pos_site_dict.keys() and infile.name.endswith("log"):
            child_logger.info("open log file {}".format(os.path.join(seq_log_folder, infile.name)))
            with open(os.path.join(seq_log_folder, infile.name), "r") as f:
                for line in f:
                    if re.search(r"-\s[{}]+$".format(required_species), line):
                        pattern = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)")
                        gene_name = (re.search(pattern, line)).group(1)
                        protein_name = (re.search(pattern, line)).group(2)
                        if not dict_of_gene_names.get(gene_name):
                            dict_of_gene_names[gene_name] = [protein_name, file_number]
    if not pattern:
        child_logger.info("log pattern not found to get gene names")
    return dict_of_gene_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?')
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?', default='')
    parser.add_argument('--required', help='Number of required (single target) species for analysis', nargs='?')
    args = parser.parse_args()
    in_dir = args.i
    log_folder = args.log
    required_species = args.required
    print("Passed args: input directory {}, log folder {}, required species {}".format(in_dir, log_folder,
                                                                                       required_species))
    try:
        main(in_dir, log_folder, required_species)
    except BaseException as e:
        print("Unexpected error: {}, \ntraceback: P{}".format(e.args, traceback.print_tb(e.__traceback__)))
    print("Common dict of genes under positive of length ", len(common_pos_gene_dict), ":\n", common_pos_gene_dict)
    print("The work has been completed")
