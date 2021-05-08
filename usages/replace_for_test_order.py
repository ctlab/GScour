#!/usr/bin/env python
import logging
import os
import argparse
import shutil


def parse_input_dir(in_folder):
    for species_folder in os.scandir(in_folder):
        if os.path.isdir(species_folder):
            for item_folder in os.scandir(species_folder):
                if os.path.isdir(item_folder):
                    yield species_folder.name, item_folder.name
                break


def main(in_folder, out_folder):
    for species_folder, item_folder in parse_input_dir(in_folder):
        shutil.copytree(os.path.join(in_folder, species_folder, item_folder),
                        os.path.join(out_folder, species_folder, item_folder))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Path to the folder with fasta files sorted by separated folders', nargs='?')
    parser.add_argument('--o', help='Path to the output folder with species folders which contain only one fasta file '
                                    'for guessing order', nargs='?')
    args = parser.parse_args()
    out_dir = args.o
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    try:
        main(args.i, out_dir)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
    print("The work has been completed")