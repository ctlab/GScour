#!/usr/bin/env python
import argparse
from scipy import stats
import os
import logging
import re
import pandas as pd

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
pattern_proportion = re.compile(r"proportion\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")  # group(1) = 0
pattern_background = re.compile(r"background\sw\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")  # group(2) = 1
pattern_foreground = re.compile(r"foreground\sw\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")  # group(3) = 2a
common_pos_gene_dict = dict()


def get_statistics(infile):
    with open(infile, "r") as f:
        pos_sites = list()
        np = 0
        ln = 0.
        background_2a, background_2b, foreground_2a, foreground_2b = "", "", "", ""
        background_0, background_1, foreground_0, foreground_1 = "", "", "", ""
        proportion_0, proportion_1, proportion_2a, proportion_2b = "", "", "", ""
        for line in f:
            if re.search(ln_np_pattern, line):
                np = int(re.search(ln_np_pattern, line).group(1))
                ln = float(re.search(ln_np_pattern, line).group(2))
            if re.search(pattern_proportion, line):
                proportion_0 = re.search(pattern_proportion, line).group(1)
                proportion_1 = re.search(pattern_proportion, line).group(2)
                proportion_2a = re.search(pattern_proportion, line).group(3)
                proportion_2b = re.search(pattern_proportion, line).group(4)
            if re.search(pattern_background, line):
                background_0 = re.search(pattern_background, line).group(1)
                background_1 = re.search(pattern_background, line).group(2)
                background_2a = re.search(pattern_background, line).group(3)
                background_2b = re.search(pattern_background, line).group(4)
            if re.search(pattern_foreground, line):
                foreground_0 = re.search(pattern_foreground, line).group(1)
                foreground_1 = re.search(pattern_foreground, line).group(2)
                foreground_2a = re.search(pattern_foreground, line).group(3)
                foreground_2b = re.search(pattern_foreground, line).group(4)
            if re.search(pos_sites_string, line):
                pattern = re.search(position_acid_pattern, f.readline())
                while pattern:
                    position = int(pattern.group(1))
                    acid = pattern.group(2)
                    probability = float(pattern.group(3))
                    pos_sites.append([position, acid, probability])
                    pattern = re.search(position_acid_pattern, f.readline())
    return np, ln, (pos_sites, proportion_0, proportion_1, proportion_2a, proportion_2b,
                    background_0, background_1, background_2a, background_2b,
                    foreground_0, foreground_1, foreground_2a, foreground_2b)


def calc_p_value(np0, ln0, np1, ln1):
    delta_ll_twice = 2 * (float(ln1) - (float(ln0)))  # twice the difference (ll from positive vs ll from neutral model)
    n = np1 - np0  # degree of freedom - difference in the number of parameters
    p_val = 1 - stats.chi2.cdf(delta_ll_twice, n)  # or stats.chisqprob(delta_ll_twice, n)
    print("n=", n, "chi2=", stats.chi2.cdf(delta_ll_twice, n), "p_val=", p_val)
    return p_val


def get_gene_name_protein_id(item_folder_name, seq_log_folder, target_species, child_logger):
    """ get gene name and protein_id which corresponds to the required_species """
    pattern, gene_name, protein_id = "", "", ""
    for infile in os.scandir(seq_log_folder):
        file_number = infile.name.split('.')[0]
        if file_number == item_folder_name and infile.name.endswith("log"):
            child_logger.info("open log file {}".format(os.path.join(seq_log_folder, infile.name)))
            with open(os.path.join(seq_log_folder, infile.name), "r") as f:
                for line in f:
                    if re.search(r"-\s[{}]+$".format(target_species), line):
                        pattern = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)")
                        try:
                            gene_name = (re.search(pattern, line)).group(1)
                            protein_id = (re.search(pattern, line)).group(2)
                        except AttributeError:
                            if not gene_name:
                                child_logger.error("Can't parse gene name from log")
                                gene_name = "gene-{}".format(file_number)
                            if not gene_name:
                                child_logger.error("Can't parse protein id from log")
                                gene_name = "protein_id-{}".format(file_number)
    if not pattern:
        child_logger.info("log pattern not found to get gene names")
    return gene_name, protein_id


def count_sites(in_folder, species_folder, item, child_logger, ortho_logs, target_species, species_folder_sheet,
                required_p_value):
    broken_paml_outs_item = list()
    no_significance_item = 0
    positive_sites_number_item = 0
    gene_protein_dict = dict()
    np0, ln0, np1, ln1 = 0, 0., 0, 0.
    pos_sites = []
    if os.path.isdir(item):
        item_folder_name = item.name
        item_id = "{}/{}".format(species_folder.name, item_folder_name)
        for infile in os.listdir(item):
            if infile.endswith("_null1_masked.out"):
                np0, ln0, _, = get_statistics(os.path.join(in_folder, species_folder.name,
                                                           item_folder_name, infile))
            if infile.endswith("_alter1_masked.out"):
                np1, ln1, sites_tuple = get_statistics(os.path.join(in_folder, species_folder.name,
                                                                    item_folder_name, infile))
                pos_sites, proportion_0, proportion_1, proportion_2a, proportion_2b, \
                background_0, background_1, background_2a, background_2b, \
                foreground_0, foreground_1, foreground_2a, foreground_2b = sites_tuple
        gene_name, protein_id = get_gene_name_protein_id(item_folder_name, ortho_logs, target_species, child_logger)
        if all([np0, np1, ln0, ln1]):
            p_val = calc_p_value(np0, ln0, np1, ln1)
            species_folder_sheet['NCBI protein_id'].append(protein_id)
            species_folder_sheet['Gene name'].append(gene_name)
            species_folder_sheet['0, proportion'].append(proportion_0)
            species_folder_sheet['0, Dn/Ds foreground'].append(foreground_0)
            species_folder_sheet['0, Dn/Ds background'].append(background_0)
            species_folder_sheet['1, proportion'].append(proportion_1)
            species_folder_sheet['1, Dn/Ds foreground'].append(foreground_1)
            species_folder_sheet['1, Dn/Ds background'].append(background_1)
            species_folder_sheet['2a, proportion'].append(proportion_2a)
            species_folder_sheet['2a, Dn/Ds foreground'].append(foreground_2a)
            species_folder_sheet['2a, Dn/Ds background'].append(background_2a)
            species_folder_sheet['2b, proportion'].append(proportion_2b)
            species_folder_sheet['2b, Dn/Ds foreground'].append(foreground_2a)
            species_folder_sheet['2b, Dn/Ds background'].append(background_2a)
            species_folder_sheet['P-value'].append(p_val)
            "the table is a main source of information, p-value and hence positive_sites_number_item, " \
            "gene_protein_dict, no_significance_item can be customized"
            if p_val and p_val < required_p_value:
                number_pos = len(pos_sites)
                child_logger.info("P.S: Item {} Gene_name {} Protein_id {} | Dn/Ds foreground (2a)={} | Dn/Ds "
                                  "foreground (2b)={}\n\t\t\t\tDn/Ds background (2a)={} |"
                                  " Dn/Ds background (2b)={}\nP-value={}, number of positive sites={}".
                                  format(item_id, gene_name, protein_id, foreground_2a, foreground_2b,
                                         background_2a, background_2b, p_val, number_pos))
                positive_sites_number_item += number_pos
                for sites in pos_sites:
                    pos, acid, probability = [sites[i] for i in range(3)]
                    child_logger.info("{} Gene_name {} positive sites : position, acid, probability : {}, {}, {}".
                                      format(item_id, gene_name, pos, acid, probability))

                if not gene_protein_dict.get(gene_name):
                    gene_protein_dict[gene_name] = (protein_id, p_val, species_folder.name)
            else:
                child_logger.info("Item {} Gene_name {} Protein_id {} | Dn/Ds foreground (2a)={} | Dn/Ds "
                                  "foreground (2b)={}\n\t\t\t\tDn/Ds background (2a)={} |"
                                  " Dn/Ds background (2b)={}\nP-value={}".
                                  format(item_id, gene_name, protein_id, foreground_2a, foreground_2b,
                                         background_2a, background_2b, p_val))
                no_significance_item += 1
        else:
            child_logger.warning("Item {} lack of params: np0 {}, ln0 {}, np1 {}, ln1 {}, pos_sites {}".format(
                item_id, np0, ln0, np1, ln1, pos_sites))
            broken_paml_outs_item.append(item_folder_name)
    return broken_paml_outs_item, no_significance_item, positive_sites_number_item, gene_protein_dict


def choose_the_lowest_p_value(joint_dict, sub_dict):
    keys_for_adding = []
    if joint_dict.keys():
        for j_key in joint_dict.keys():
            for s_key in sub_dict.keys():
                if j_key == s_key:
                    if sub_dict[s_key][1] < joint_dict[s_key][1]:
                        joint_dict[s_key] = (sub_dict[s_key][0], sub_dict[s_key][1],
                                             "{},{}".format(sub_dict[s_key][2], joint_dict[s_key][2]))
                    else:
                        joint_dict[j_key] = (joint_dict[j_key][0], joint_dict[j_key][1],
                                             "{},{}".format(joint_dict[j_key][2], sub_dict[j_key][2]))
                else:
                    keys_for_adding.append(s_key)
        for key in keys_for_adding:
            joint_dict[key] = (sub_dict[key][0], sub_dict[key][1], sub_dict[key][2])
    else:
        joint_dict.update(sub_dict)


def main(in_folder, ortho_logs, target_species, required_p_value):
    global common_pos_gene_dict
    for species_folder in os.scandir(in_folder):
        species_folder_sheet = {
            'NCBI protein_id': [], 'Gene name': [], '0, proportion': [],
            '0, Dn/Ds foreground': [], '0, Dn/Ds background': [],
            '1, proportion': [], '1, Dn/Ds foreground': [], '1, Dn/Ds background': [],
            '2a, proportion': [], '2a, Dn/Ds foreground': [], '2a, Dn/Ds background': [],
            '2b, proportion': [], '2b, Dn/Ds foreground': [], '2b, Dn/Ds background': [],
            'P-value': []
            }
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
                                                    ortho_logs, target_species, species_folder_sheet, required_p_value)

                    broken_paml_outs += broken_paml_outs_item
                    no_significance += no_significance_item
                    positive_sites_number += positive_sites_number_item

                    choose_the_lowest_p_value(genes_under_positive, gene_protein_dict)
                    print("genes_under_positive for {} of len {}, {}".format(species_folder.name,
                                                                             len(genes_under_positive),
                                                                             genes_under_positive))
                    # genes_under_positive.update(gene_protein_dict)

            child_logger.warning("Species folder {}: broken_paml_outs : {} : {}".
                                 format(species_folder.name, len(broken_paml_outs), broken_paml_outs))
            child_logger.info("Species folder {}: number of no significance files {}".format(species_folder.name,
                                                                                             no_significance))
            child_logger.info("Species folder {}: number of positive sites {}".format(species_folder.name,
                                                                                      positive_sites_number))
            child_logger.info(
                "Species folder {} Number of genes under P.S={}: gene_name : protein_id \n{}".format(
                    species_folder.name, len(genes_under_positive), repr(genes_under_positive)))

            choose_the_lowest_p_value(common_pos_gene_dict, genes_under_positive)
            print("common_pos_gene_dict", common_pos_gene_dict)
            # common_pos_gene_dict.update(genes_under_positive)
            df = pd.DataFrame(species_folder_sheet, columns=['NCBI protein_id', 'Gene name', '0, proportion',
                                                             '0, Dn/Ds foreground', '0, Dn/Ds background',
                                                             '1, proportion', '1, Dn/Ds foreground',
                                                             '1, Dn/Ds background',
                                                             '2a, proportion', '2a, Dn/Ds foreground',
                                                             '2a, Dn/Ds background',
                                                             '2b, proportion', '2b, Dn/Ds foreground',
                                                             '2b, Dn/Ds background',
                                                             'P-value'])

            df.to_excel(writer, sheet_name=species_folder.name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?', required=True)
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?', required=True)
    parser.add_argument('--required', help='Number of required (single target) species for analysis', nargs='?',
                        required=True)
    parser.add_argument('--p', help='p-value level', nargs='?', required=True)
    args = parser.parse_args()
    in_dir = args.i
    log_folder = args.log
    required_species = args.required
    p_value_required = float(args.p)
    print("Passed args: input directory {}, log folder {}, required species {}, p-value {}".format(in_dir, log_folder,
                                                                                                   required_species,
                                                                                                   p_value_required))

    try:
        common_sheet_path = os.path.join(in_dir, 'masked_common_sheet.xlsx')
        writer = pd.ExcelWriter(common_sheet_path, engine='xlsxwriter')
        main(in_dir, log_folder, required_species, p_value_required)
        values = list(common_pos_gene_dict.values())
        # print("value[0],\n", values[0], "value[1],\n", values[1],
        #       "\nentire keys\n", common_pos_gene_dict.keys())
        print("Results are recorded in {}".format(common_sheet_path))
        summary_sheet = {
            'Gene name': list(common_pos_gene_dict.keys()), 'NCBI protein_id':
                [i[0] for i in values], 'p-value': [i[1] for i in values], 'Species group': [i[2] for i in values]
            }
        df = pd.DataFrame(summary_sheet, columns=['Gene name', 'NCBI protein_id', 'p-value', 'Species group'])
        df.to_excel(writer, sheet_name='summary')
        writer.save()
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
    print("Common dict of genes under positive of length ", len(common_pos_gene_dict), ":\n", common_pos_gene_dict)
    print("The work has been completed")
