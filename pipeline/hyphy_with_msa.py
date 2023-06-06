#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import json
import traceback
from subprocess import Popen, PIPE
import os
import re
import pandas as pd

result = {'FEL_positive': list(), 'FEL_negative': list(), 'BUSTED_positive': list(), 'ABSREL_positive': list(),
          'MEME_positive': list(), 'RELAX': list(), 'RELAX_selection_intensity_parameter': list(),
          'RELAX_selection_strength': list()}
positive = 0
negative = 0
relax = 0
exception_analysis = set()
exception_parsing = set()
analyzed_items = 0


def parse_fel(gene_folder_name, out_file):
    global positive
    global negative
    """### ** Found _0_ sites under pervasive positive diversifying and _8_ sites under negative selection at p <=
    0.05** """
    try:
        with open(out_file, 'r') as f:
            last_line = f.readlines()[-1]
            print('last line', last_line)
            try:
                if 'Check errors.log for execution error details.' in last_line:
                    exception_analysis.add(gene_folder_name)
                    return
                pos = int(re.search('###\s\*\*\sFound\s_(\d+)_\ssites\sunder\spervasive\spositive\sdiversifying', last_line).group(1))
                neg = int(re.search('and\s_(\d+)_\ssites\sunder\snegative\sselection', last_line).group(1))
                if pos:
                    result['FEL_positive'].append(gene_folder_name.split('_')[0])
                    positive += 1
                if neg:
                    result['FEL_negative'].append(gene_folder_name.split('_')[0])
                    negative += 1
            except AttributeError:
                print(f'error FEL, not such string in {gene_folder_name}/{out_file}')
                exception_parsing.add(gene_folder_name)
    except Exception:
        print(traceback.format_exc())


def launch_fel(gene_folder_name, al, tree, relaunch):
    output_file = 'FEL.out'
    if relaunch == 'y' or (relaunch == 'n' and not output_file):
        with open(output_file, 'w') as file:
            print('1st al FEL', al)
            proc = Popen([f'hyphy pre-msa.bf --input {al}'], shell=True, stdout=file,
                         stderr=file)
            proc.wait()
            # print('alignment_base_name', alignment_base_name)
            proc = Popen([f'muscle -align {al}_protein.fas -output {al}_protein.msa'], shell=True, stdout=file,
                         stderr=file)
            proc.wait()
            proc = Popen(
                [f'hyphy post-msa.bf --protein-msa {al}_protein.msa --nucleotide-sequences '
                 f'{al}_nuc.fas --output {al}.msa --compress No'],
                shell=True, stdout=file,
                stderr=file)
            proc.wait()
            proc = Popen([f'hyphy fel --alignment {al}.msa --tree {tree} --pvalue 0.05'], shell=True, stdout=file,
                         stderr=file)
            proc.wait()
        if proc.returncode != 0:
            print(f'An error occurred during FEL execution with {gene_folder_name}')
            exception_analysis.add(gene_folder_name)
            return

    print(f"Start parsing output file '{output_file}'.")
    parse_fel(gene_folder_name, output_file)


def parse_meme(gene_folder_name, out_file):
    """ ### ** Found _0_ sites under episodic diversifying positive selection at p <= 0.05** """
    global positive
    try:
        with open(out_file, 'r') as f:
            last_line = f.readlines()[-1]
            print('last line', last_line)
            try:
                if 'Check errors.log for execution error details.' in last_line:
                    exception_analysis.add(gene_folder_name)
                    return
                pos = int(re.search('###\s\*\*\sFound\s_(\d+)_\ssites\sunder\sepisodic\sdiversifying\spositive', last_line).group(1))
                if pos != 0:
                    result['MEME_positive'].append(gene_folder_name.split('_')[0])
                    positive += 1
            except AttributeError:
                print(f'error MEME, not such string in {gene_folder_name}/{out_file}')
                exception_parsing.add(gene_folder_name)
    except Exception:
        print(traceback.format_exc())
        exception_parsing.add(gene_folder_name)


def launch_meme(gene_folder, al, tree, relaunch):
    output_file = 'MEME.out'
    if relaunch == 'y' or (relaunch == 'n' and not output_file):
        with open(output_file, 'w') as file:
            proc = Popen([f'hyphy meme --alignment {al} --tree {tree} --pvalue 0.05'], shell=True, stdout=file,
                         stderr=file)
            proc.wait()
        if proc.returncode != 0:
            print(f'An error occurred during MEME execution with {gene_folder}')
            exception_analysis.add(gene_folder)
            return
    print(f"Start parsing output file '{output_file}'.")
    parse_meme(gene_folder, output_file)


def parse_absrel(gene_folder_name, al):
    global positive
    json_path = f'{al}.ABSREL.json'
    print('alignment name ', al)
    print(f"Start parsing output file {json_path}")
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
            if data['test results']['positive test results'] != 0:
                result['ABSREL_positive'].append(gene_folder_name.split('_')[0])
                positive += 1
                print('test results', data['test results'])
    except Exception:
        print(traceback.format_exc())
        exception_parsing.add(gene_folder_name)


def launch_absrel(gene_folder, al, tree, relaunch):
    output_file = 'ABSREL.out'
    if relaunch == 'y' or (relaunch == 'n' and not output_file):
        with open(output_file, 'w') as file:
            proc = Popen([f'hyphy absrel --alignment {al} --tree {tree} --pvalue 0.05'], shell=True, stdout=file,
                         stderr=file)
            proc.wait()
        if proc.returncode != 0:
            print(f'An error occurred during ABSREL execution in {gene_folder}')
            exception_analysis.add(gene_folder)
            return
    parse_absrel(gene_folder, al)


def parse_busted(gene_folder_name, al):
    global positive
    json_path = f'{al}.BUSTED.json'
    print('align name ', al)
    print(f"Start parsing output file '{json_path}'.")
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
            if data['test results']['p-value'] < 0.05:
                result['BUSTED_positive'].append(gene_folder_name.split('_')[0])
                positive += 1
                print('test results', data['test results'])
    except Exception:
        print(traceback.format_exc())
        exception_parsing.add(gene_folder_name)


def launch_busted(gene_folder, al, tree, relaunch):
    output_file = 'BUSTED.out'
    if relaunch == 'y' or (relaunch == 'n' and not output_file):
        with open(output_file, 'w') as file:
            proc = Popen([f'hyphy busted --alignment {al} --tree {tree} --pvalue 0.05'], shell=True, stdout=file,
                         stderr=file)
            proc.wait()
        if proc.returncode != 0:
            print(f'An error occurred during BUSTED execution with {gene_folder}')
            exception_analysis.add(gene_folder)
            return
    parse_busted(gene_folder, al)


def parse_relax(gene_folder_name, al):
    global relax
    json_path = f'{al}.RELAX.json'
    print(f"Start parsing output file '{json_path}'.")
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
            if data['test results']['p-value'] < 0.05:
                result['RELAX'].append(gene_folder_name.split('_')[0])
                intensity_parameter = data['test results']['relaxation or intensification parameter']
                result['RELAX_selection_intensity_parameter'].append(intensity_parameter)
                if intensity_parameter > 1:
                    result['RELAX_selection_strength'].append('intensified')
                if intensity_parameter < 1:
                    result['RELAX_selection_strength'].append('relaxed')
                relax += 1
                print('test results', data['test results'])
    except Exception:
        print(traceback.format_exc())
        exception_parsing.add(gene_folder_name)


def launch_relax(gene_folder, al, tree, relaunch):
    output_file = 'RELAX.out'
    if relaunch == 'y' or (relaunch == 'n' and not output_file):
        with open(output_file, 'w') as file:
            proc = Popen([
                f'hyphy relax --alignment {al} --tree {tree} --pvalue 0.05]'],
                shell=True, stdin=PIPE, stdout=file, stderr=PIPE, bufsize=1,
                universal_newlines=True)
            proc.communicate(input="2\n")  # for foreground - human
            proc.wait()
        if proc.returncode != 0:
            print(f'An error occurred during RELAX execution with {gene_folder}')
            exception_analysis.add(gene_folder)
            return
    parse_relax(gene_folder, al)


def main(input_folder, extension, relaunch):
    global analyzed_items
    for gene_folder in os.scandir(input_folder):
        alignment, tree, tag_tree = '', '', ''
        if os.path.isdir(gene_folder):
            os.chdir(gene_folder)
            print(f'Working with {gene_folder.name}')
            for f in os.scandir(gene_folder):
                if f.name.split(f'{gene_folder.name}.')[-1] == extension:
                    alignment = os.path.join(input_folder, gene_folder.name, f.name)
                if f.name.split('_')[-1] == 'tag.tree':  # great_apes_tag.tree
                    tag_tree = os.path.join(input_folder, gene_folder.name, f.name)
                if f.name.split('_')[-1] == 'apes.tree':
                    tree = os.path.join(input_folder, gene_folder.name, f.name)
            if alignment and tree and tag_tree:
                launch_fel(gene_folder.name, alignment, tree, relaunch)
                msa_alignment = f'{alignment}.msa'
                launch_meme(gene_folder.name, msa_alignment, tree, relaunch)
                launch_absrel(gene_folder.name, msa_alignment, tree, relaunch)
                launch_busted(gene_folder.name, msa_alignment, tree, relaunch)
                launch_relax(gene_folder.name, msa_alignment, tag_tree, relaunch)
            else:
                print('No necessary files')
        analyzed_items += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i',
                        help='The full path to the folder contains items(genes) folders with input files for hyphy',
                        nargs='?', required=True)
    parser.add_argument('--e',
                        help='Extension for input files for',
                        nargs='?', default='best.nuc.fas')
    parser.add_argument('--r',
                        help='Rework: \'y\'-relaunch alignments and methods, \'n\'-gather results where they are and '
                             'launch methods and gather where not',
                        nargs='?', default='n')
    args = parser.parse_args()
    try:
        main(args.i, args.e, args.r)
    except Exception:
        print(traceback.format_exc())
    finally:
        df = pd.DataFrame.from_dict(result, orient='index')
        df = df.transpose()
        print(df)
        print(result)
        os.chdir(args.i)
        df.to_csv('hyphy_results.csv')
        print(f'Summary: positive = {positive}, negative = {negative}, relax={relax}, '
              f'exception_analysis = {len(exception_analysis)}: {exception_analysis},'
              f' exception_parsing = {len(exception_parsing)}: {exception_parsing}, analyzed_items = {analyzed_items}')
        print(f'Results are in {args.i}/hyphy_results.csv. Done')
