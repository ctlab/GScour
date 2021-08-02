#!/usr/bin/env python
import argparse
import os
import re
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def get_tree_and_codes(rst_file_path, trees_folder, species_folder_name):
    tree_str, codes_str = "", ""
    # tree_node_labels_str = "", "", ""
    with open(os.path.join(trees_folder, "{}.tree".format(species_folder_name)), 'r') as tf:
        tree_str = tf.readline()
    with open(rst_file_path, 'r') as f:
        counter = 0
        # counter_node_labels = 0
        for line in f:
            # if counter_node_labels == 1:
            #     tree_node_labels_str = line
            #     break
            if counter == 1:
                counter += 1
            elif counter == 2:
                codes_str = line
                break
                # counter += 1
                # continue
            if re.search('^\([^:]+\);$', line):
                # tree_str = line
                counter += 1
            # if re.match('tree with node labels', line):
            #     counter_node_labels += 1
    return tree_str, codes_str  # , tree_node_labels_str


def collect_previous_tree_items(tree_string, i):
    tree_items_list = []
    for tree_item in tree_string.split()[i::-1]:
        numeric_species_name = re.search(r'\d+', tree_item).group(0)
        tree_items_list.append(numeric_species_name)
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


# def replace_for_true_labels(branchcodes, node_labels):
#     logger.info("branchcodes {}".format(repr(branchcodes)))
#     logger.info("node_labels {}".format(repr(node_labels)))
#     labels_dict = {}
#     for line in branchcodes.split('\n'):
#         logger.info("line {}".format(line))
#         if re.search(r'\s(\d+)$', line):
#             label = re.search(r'\s(\d+)$', line).group(1)
#             logger.info("label {}".format(repr(label)))
#             true_label = re.search(label + r'\_(\d+)', node_labels).group(1)
#             labels_dict[label] = true_label
#             logger.info("true_label {}".format(true_label))
#             edited_line = re.sub(r'(\d+)$', true_label, line)
#             logger.info("edited line {}".format(repr(edited_line)))
#             branchcodes = branchcodes.replace(line, edited_line)
#             logger.info("branchcodes replaced {}".format(repr(branchcodes)))
#         if re.search(r'\s((\d,)+\d+)$', line):
#             labels = re.search(r'\s((\d,)+\d+)$', line).group(1)
#             true_labels = labels
#             i = 0
#             for label in labels.split(','):
#                 if i != len(labels.split(',')):
#                     true_labels = re.sub(label + ',', labels_dict[label] + ',', true_labels)
#                     i += 1
#                 else:
#                     true_labels = re.sub(label, labels_dict[label], true_labels)
#                 logger.info("multi: label {}, true labels {}".format(repr(label), true_labels))
#             edited_line = re.sub(labels, true_labels, line)
#             logger.info("edited line {}".format(edited_line))
#             branchcodes = branchcodes.replace(line, edited_line)
#             logger.info("branchcodes replaced {}".format(branchcodes))
#     logger.info("edited branchcodes {}".format(branchcodes))
#     return branchcodes


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
                        # phylogenetic_tree, codes_for_branches, node_labels = get_tree_and_codes(rst_file_path)
                        phylogenetic_tree, codes_for_branches = get_tree_and_codes(rst_file_path, trees_folder,
                                                                                   species_folder.name)
                        logger.info("phylogenetic_tree {}\ncodes_for_branches {}".format(phylogenetic_tree,
                                                                                         codes_for_branches))
                        if phylogenetic_tree and codes_for_branches:
                            branchcodes = get_branchcodes(phylogenetic_tree, codes_for_branches)
                            # branchcodes_true_labels = replace_for_true_labels(branchcodes, node_labels)
                            # logger.info("AFTER replacing branchcodes_true_labels {}, node_labels {}".
                            #             format(branchcodes_true_labels, node_labels))
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
