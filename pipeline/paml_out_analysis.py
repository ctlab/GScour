#!/usr/bin/env python
import argparse
import sys
from scipy import stats
import os
import logging
import re

"""
Site class 2a: Codon sites evolving under positive selection in the selected branch (dN/dS>1),
and under purifying selection in the rest of the tree (0<dN/dS<1)
Site class 2b: Codon sites evolving under positive selection in the selected branch (dN/dS>1),
and under neutral evolution in the rest of the tree (dN/dS =1) 
"""
ln_np_pattern = re.compile(r"lnL\(ntime:\s+\d+\s+np:\s+(\d+)\):\s+(-\d+\.\d+)")
# POSITIVE_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
pos_sites_string = re.compile(r"Positive\ssites\sfor\sforeground\slineages\sProb\(w>1\):")
position_acid_pattern = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
pattern_background = re.compile(r"background\sw\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)")  # group(1) = 2a
pattern_foreground = re.compile(r"foreground\sw\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)")  # group(2) = 2b
common_pos_gene_dict = dict()


def get_statistics(infile):
    with open(infile, "r") as f:
        pos_sites = list()
        np = 0
        ln = 0.
        for line in f:
            if re.search(ln_np_pattern, line):
                np = int(re.search(ln_np_pattern, line).group(1))
                ln = float(re.search(ln_np_pattern, line).group(2))
            if re.search(pattern_background, line):
                background_2a = re.search(pattern_background, line).group(1)
                background_2b = re.search(pattern_background, line).group(2)
            if re.search(pattern_foreground, line):
                foreground_2a = re.search(pattern_foreground, line).group(1)
                foreground_2b = re.search(pattern_foreground, line).group(2)
            if re.search(pos_sites_string, line):
                pattern = re.search(position_acid_pattern, f.readline())
                while pattern:
                    position = int(pattern.group(1))
                    acid = pattern.group(2)
                    probability = float(pattern.group(3))
                    pos_sites.append([position, acid, probability])
                    pattern = re.search(position_acid_pattern, f.readline())
    return np, ln, pos_sites, background_2a, background_2b, foreground_2a, foreground_2b


def calc_p_value(np0, ln0, np1, ln1):
    ln_val = 2 * (float(ln1) - (float(ln0)))
    n = np1 - np0
    p_val = 1 - stats.chi2.cdf(ln_val, n)
    return p_val


def get_gene_name_protein_id(item_folder_name, seq_log_folder, target_species, child_logger):
    """ get gene name and protein_id wich corresponds to the required_species """
    pattern, gene_name, protein_id = "", "", ""
    for infile in os.scandir(seq_log_folder):
        file_number = infile.name.split('.')[0]
        if file_number == item_folder_name and infile.name.endswith("log"):
            child_logger.info("open log file {}".format(os.path.join(seq_log_folder, infile.name)))
            with open(os.path.join(seq_log_folder, infile.name), "r") as f:
                for line in f:
                    if re.search(r"-\s[{}]+$".format(target_species), line):
                        pattern = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)")
                        gene_name = (re.search(pattern, line)).group(1)
                        protein_id = (re.search(pattern, line)).group(2)
    if not pattern:
        child_logger.info("log pattern not found to get gene names")
    return gene_name, protein_id


def count_sites(in_folder, species_folder, item, child_logger, ortho_logs, target_species):
    broken_paml_outs_item = list()
    no_significance_item = 0
    positive_sites_number_item = 0
    gene_protein_dict_item = dict()

    np0, ln0, np1, ln1 = 0, 0., 0, 0.
    background_2a, background_2b, foreground_2a, foreground_2b = "", "", "", ""
    pos_sites = []
    if os.path.isdir(item):
        item_folder_name = item.name
        item_id = "{}/{}".format(species_folder.name, item_folder_name)
        for infile in os.listdir(item):
            if infile.endswith("_null1_masked.out"):
                np0, ln0, _, _, _, _, _ = get_statistics(os.path.join(in_folder, species_folder.name,
                                                                      item_folder_name, infile))
            if infile.endswith("_alter1_masked.out"):
                np1, ln1, pos_sites, background_2a, background_2b, foreground_2a, foreground_2b \
                    = get_statistics(os.path.join(in_folder, species_folder.name,
                                                  item_folder_name, infile))

        gene_name, protein_id = get_gene_name_protein_id(item_folder_name, ortho_logs, target_species, child_logger)
        if all([np0, np1, ln0, ln1]):
            p_val = calc_p_value(np0, ln0, np1, ln1)
            if p_val and p_val < 0.05:
                number_pos = len(pos_sites)
                child_logger.info("P.S: Item {} Gene_name {} Protein_id {} | Dn/Ds Target species (2a)={} | Dn/Ds "
                                  "Target species (2b)={}\n\t\t\t\tDn/Ds background species (2a)={} |"
                                  " Dn/Ds background species (2b)={}\nP-value={}, number of positive sites={}".
                                  format(item_id, gene_name, protein_id, foreground_2a, foreground_2b,
                                         background_2a, background_2b, p_val, number_pos))
                positive_sites_number_item += number_pos
                for sites in pos_sites:
                    pos, acid, probability = [sites[i] for i in range(3)]
                    child_logger.info("{} Gene_name {} positive sites : position, acid, probability : {}, {}, {}".
                                      format(item_id, gene_name, pos, acid, probability))

                if not gene_protein_dict_item.get(gene_name):
                    gene_protein_dict_item[gene_name] = [protein_id, item_folder_name]
            else:
                child_logger.info("Item {} Gene_name {} Protein_id {} | Dn/Ds Target species (2a)={} | Dn/Ds "
                                  "Target species (2b)={}\n\t\t\t\tDn/Ds background species (2a)={} |"
                                  " Dn/Ds background species (2b)={}\nP-value={}".
                                  format(item_id, gene_name, protein_id, foreground_2a, foreground_2b,
                                         background_2a, background_2b, p_val))
                no_significance_item += 1
        else:
            child_logger.warning("Item {} lack of params: np0 {}, ln0 {}, np1 {}, ln1 {}, pos_sites {}".format(
                item_id, np0, ln0, np1, ln1, pos_sites))
            broken_paml_outs_item.append(item_folder_name)
    return broken_paml_outs_item, no_significance_item, positive_sites_number_item, gene_protein_dict_item


def main(in_folder, ortho_logs, target_species):
    global common_pos_gene_dict
    for species_folder in os.scandir(in_folder):
        broken_paml_outs = list()
        no_significance = 0
        positive_sites_number = 0
        genes_under_positive = dict()

        if os.path.isdir(species_folder):
            result_file = os.path.join(in_folder, species_folder.name, "{}.{}".format(species_folder.name,
                                                                                      "result"))
            child_logger = logging.getLogger('__main__.' + species_folder.name)
            child_logger.addHandler(logging.FileHandler(result_file))
            child_logger.setLevel(10)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    broken_paml_outs_item, no_significance_item, positive_sites_number_item, \
                        gene_protein_dict = count_sites(in_folder, species_folder, item, child_logger,
                                                        ortho_logs, target_species)

                    broken_paml_outs += broken_paml_outs_item
                    no_significance += no_significance_item
                    positive_sites_number += positive_sites_number_item
                    genes_under_positive.update(gene_protein_dict)

            child_logger.warning("Species folder {}: broken_paml_outs : {} : {}".
                                 format(species_folder.name, len(broken_paml_outs), broken_paml_outs))
            child_logger.info("Species folder {}: number of no significance files {}".format(species_folder.name,
                                                                                             no_significance))
            child_logger.info("Species folder {}: number of positive sites {}".format(species_folder.name,
                                                                                      positive_sites_number))
            child_logger.info(
                "Species folder {} Number of genes under P.S={}: gene_name : protein_id : item \n{}".format(
                    species_folder.name, len(genes_under_positive), repr(genes_under_positive)))

            for key in genes_under_positive.keys():
                if key not in common_pos_gene_dict.keys():
                    common_pos_gene_dict[key] = list()
                common_pos_gene_dict[key] += genes_under_positive[key]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?')
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?')
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
        print("Unexpected error: {}".format(e))
    print("Common dict of genes under positive of length ", len(common_pos_gene_dict), ":\n", common_pos_gene_dict)
    print("The work has been completed")
