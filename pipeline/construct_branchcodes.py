#!/usr/bin/env python
import argparse
import os
import re
import logging

"""construct branchcodes for SWAMP"""

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def get_tree_and_codes(rst_file_path, trees_folder, species_folder_name):
    tree_str, codes_str = "", ""
    with open(os.path.join(trees_folder, "{}.tree".format(species_folder_name)), 'r') as tf:
        tree_str = tf.readline()
    with open(rst_file_path, 'r') as f:
        counter = 0
        for line in f:
            if counter == 1:
                counter += 1
            elif counter == 2:
                codes_str = line
                break
            if re.search('^\([^:]+\);$', line):
                counter += 1
    return tree_str, codes_str


def collect_previous_tree_items(tree_string, i):
    tree_items_list = []
    print("tree_string", tree_string)
    for tree_item in tree_string.split()[i::-1]:
        numeric_species_name = re.search(r'\d+', tree_item).group(0)
        tree_items_list.append(numeric_species_name)
        print("tree_items_list", tree_items_list)
        if re.search(r'\(', tree_item):
            removable_idx = str.find(tree_item, '(')
            new_tree_item = tree_item[0:removable_idx] + tree_item[removable_idx + 1:]
            tree_string_changed = tree_string.replace(tree_item, new_tree_item)
            tree_items_string = ",".join(tree_items_list)
            break
    return tree_items_string, tree_string_changed


def get_branchcodes(tree_string, codes_string):
    i = 0
    prev_tree_items = []
    branchcodes_list = []
    for code in codes_string.split():
        if prev_tree_items:
            branchcodes_list.append(' '.join([code, prev_tree_items[0]]))
            prev_tree_items.pop(0)
            continue
        if i < len(tree_string.split()):
            numeric_species_name = re.search(r'\d', tree_string.split()[i]).group(0)
            branchcodes_list.append(' '.join([code, numeric_species_name]))
            closing_brackets = [pos.start() for pos in re.finditer(r'\)', tree_string.split()[i])]
            if closing_brackets:
                for j in closing_brackets:
                    prev_tree_item, tree_string = collect_previous_tree_items(tree_string, i)
                    prev_tree_items.append(prev_tree_item)
            i += 1
        else:
            i -= 1
            prev_tree_item, tree_string = collect_previous_tree_items(tree_string, i)
            prev_tree_items.append(prev_tree_item)
            branchcodes_list.append(' '.join([code, prev_tree_items[0]]))
    return "\n".join(branchcodes_list)


def find_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def main(input_folder, trees_folder, branchcodes_folder):
    positive_result = []
    list_of_folders = []
    for species_folder in os.scandir(input_folder):
        if os.path.isdir(species_folder):
            logger.info("species folder {}".format(species_folder))
            list_of_folders.append(species_folder.name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    rst_file_path = find_file('rst', os.path.join(input_folder, species_folder.name, item.name))
                    if rst_file_path:
                        logger.info("rst_file_path found {}".format(rst_file_path))
                        phylogenetic_tree, codes_for_branches = get_tree_and_codes(rst_file_path, trees_folder,
                                                                                   species_folder.name)
                        logger.info("phylogenetic_tree {}\ncodes_for_branches {}".format(phylogenetic_tree,
                                                                                         codes_for_branches))
                        if phylogenetic_tree and codes_for_branches:
                            branchcodes = get_branchcodes(phylogenetic_tree, codes_for_branches)
                            with open(os.path.join(branchcodes_folder, "{}.{}".format(species_folder.name, "code")),
                                      "w") as f:
                                f.write(branchcodes)
                                logger.info("Branchcode for species_folder {} recorded:\n{}".
                                            format(species_folder, branchcodes))
                                positive_result.append(species_folder.name)
                                break
    lack_of_codes = list(set(list_of_folders) - set(positive_result))
    if lack_of_codes:
        logger.warning("Phylogenetic tree and codes for branches has not been found for species folders {}".
                       format(lack_of_codes))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?', required=True)
    parser.add_argument('--t', help='The full path to the folder contains trees',
                        nargs='?', required=True)
    parser.add_argument('--b', help='The full path to a folder for recording branchcodes',
                        nargs='?', required=True)

    args = parser.parse_args()
    input_dir = args.i
    branchcodes_dir = args.b
    trees_dir = args.t
    if not os.path.isdir(branchcodes_dir):
        os.makedirs(branchcodes_dir)
    logger.info("Passed args: input directory {}, trees directory {}, branchcode's directory {}".format(input_dir,
                                                                                                        trees_dir,
                                                                                                        branchcodes_dir))
    try:
        main(input_dir, trees_dir, branchcodes_dir)
    except BaseException as e:
        logger.warning("Unexpected error: {}".format(e))
        raise e
    logger.info("The work has been completed")
