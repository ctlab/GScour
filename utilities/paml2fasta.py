#!/usr/bin/env python
import argparse
import os
import re
"""
the script converts paml sequence format to fasta
"""


def get_input_phy_file(folder_in, masked):
    """ parse root folder with files for paml """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder) and species_folder.name.isdigit():
            print("working with species folder {}".format(species_folder.name))
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if masked == 'y':
                            if infile.endswith("_masked.phy"):
                                yield folder_in, species_folder.name, item.name, infile
                        if masked == 'n':
                            if infile.split('.')[-1] == 'phy' and not re.search('[a-zA-Z]', infile.split('.')[0]):
                                yield folder_in, species_folder.name, item.name, infile


def paml2fasta(in_file, out_file):
    with open(in_file, 'r') as fi:
        with open(out_file, 'a') as fo:
            i = 0
            for line in fi:
                if i == 0:
                    i += 1
                    continue
                if i % 2 == 1:
                    seq_name = line.strip()
                    i += 1
                    continue
                if i % 2 == 0:
                    seq = line.strip()
                    fo.write(">" + seq_name + "\n" + seq + "\n")
                    i += 1


def main(in_folder, out_folder, masked_flag):
    for folder_in, species_folder, item, infile in get_input_phy_file(in_folder, masked_flag):
        in_file = os.path.join(folder_in, species_folder, item, infile)
        out_file = os.path.join(out_folder, infile.replace('phy', 'fasta'))
        paml2fasta(in_file, out_file)
        print("file {} was written".format(out_file))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?', required=True)
    parser.add_argument('--o', help='Path to the folder with result fasta files, if it does not exist, '
                                    'it will be created automatically', nargs='?', required=True)
    parser.add_argument('--m', help='Masked "y" or "n" files to convert', nargs='?', required=True)
    args = parser.parse_args()
    out_dir = args.o
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    main(args.i, out_dir, args.m)
    print("The work has been completed")