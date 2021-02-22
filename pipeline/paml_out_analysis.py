#!/usr/bin/env python
import argparse
import sys
import logging
from scipy import stats
import os
import logging
import re

LN_NP_PATTERN = re.compile(r"lnL\(ntime:\s+\d+\s+np:\s+(\d+)\):\s+(-\d+\.\d+)")
POSITIVE_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
POS_SITES_STRING = re.compile(r"Positive\ssites\sfor\sforeground\slineages\sProb\(w>1\):")
POSITION_ACID_PATTERN = re.compile(r"\s*(\d+)\s+(\w)\s+(\d+\.\d+)")
LOG_FILE = "analyse_paml_out.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
BROKEN_PAML_OUTS = list()
NO_SIGNIFICANCE = 0
POSITIVE_SITES_NUMBER = 0
POSITIVE_SITES_DICT = dict()
POSITIVE_GENES = list()


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


def main(folder_in):
    global NO_SIGNIFICANCE
    global POSITIVE_SITES_NUMBER
    for personal_folder in os.scandir(folder_in):
        if os.path.isdir(personal_folder):
            folder_name = personal_folder.name
            np0, ln0, np1, ln1 = 0, 0., 0, 0.
            for infile in os.listdir(personal_folder):
                if infile.endswith("_null1.out"):
                    np0, ln0, _ = get_ln_np(os.path.join(folder_in, personal_folder, infile))
                if infile.endswith("_alter1.out"):
                    np1, ln1, pos_sites = get_ln_np(os.path.join(folder_in, personal_folder, infile))
            if all([np0, np1, ln0, ln1]):
                p_val = calc_p_value(np0, ln0, np1, ln1)
                if p_val and p_val < 0.05:
                    number_pos = len(pos_sites)
                    logging.info("{} number of positive sites {}".format(folder_name, number_pos))
                    POSITIVE_SITES_NUMBER += number_pos
                    for sites in pos_sites:
                        pos, acid, probability = [sites[i] for i in range(3)]
                        logging.info("{} positive sites : position, acid, probability : {}, {}, {}".format(folder_name,
                                                                                                           pos,
                                                                                                           acid,
                                                                                                           probability))
                        if folder_name not in POSITIVE_GENES:
                            POSITIVE_GENES.append(folder_name)
                        if not POSITIVE_SITES_DICT.get(folder_name):
                            POSITIVE_SITES_DICT[folder_name] = pos_sites
                else:
                    logging.info("{} no significance, p-value {}".format(folder_name, p_val))
                    NO_SIGNIFICANCE += 1
            else:
                logging.warning("lack of params {}: np0 {}, ln0 {}, np1 {}, ln1 {}, pos_sites {}".format(
                    folder_name, np0, ln0, np1, ln1, pos_sites))
                BROKEN_PAML_OUTS.append(folder_name)


def get_genes_under_positive(genes_under_positive, log_folder):
    dict_of_gene_names = dict()
    for infile in os.listdir(log_folder):
        file_number = infile.split('.')[0]
        if file_number in genes_under_positive and infile.endswith("log"):
            with open(os.path.join(log_folder, infile), "r") as f:
                for line in f:
                    if re.search(r"-\s5$", line):
                        pattern = re.compile(r"^([a-zA-Z0-9]+)\s-\s([a-zA-Z0-9_\.]+)")
                        gene_name = (re.search(pattern, line)).group(1)
                        protein_name = (re.search(pattern, line)).group(2)
                        if not dict_of_gene_names.get(infile):
                            dict_of_gene_names[file_number] = [gene_name, protein_name]
    return dict_of_gene_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='Path to the folder with input files for paml', nargs='?')
    parser.add_argument('--log', help='Path to the log folder of "get_ortho_nucleotides.py"', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder)
        if BROKEN_PAML_OUTS:
            logging.warning("BROKEN_PAML_OUTS : {} : {}".format(len(BROKEN_PAML_OUTS), BROKEN_PAML_OUTS))
        logging.info("Number of no significance files {}".format(NO_SIGNIFICANCE))
        logging.info("Number of positive sites {}".format(POSITIVE_SITES_NUMBER))
        logging.info("Number of positive genes {} : {}".format(len(POSITIVE_GENES), POSITIVE_GENES))
        logging.info("Positive sites : file : position, acid, probability\n{}".format(repr(POSITIVE_SITES_DICT)))
        gene_names_dict = get_genes_under_positive(POSITIVE_GENES, args.log)
        if gene_names_dict:
            logging.info("Genes under positive selection: file: gene, protein:\n{}".format(repr(gene_names_dict)))
    except:
        logging.exception("Unexpected error")
    logging.info("The work has been completed")