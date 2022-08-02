#!/usr/bin/env python
import argparse
import subprocess
import sys
import functools, multiprocessing, inspect
import pandas as pd
from Bio.Phylo.PAML import codeml
import logging
import os

"""
Under this model, the relationship holds that ω = dN/dS, the ratio of
nonsynonymous/synonymous substitution rates. This basic model is fitted by specifying model = 0
NSsites = 0
The ω ratio is a measure of natural selection acting on the protein. Simplistically, values of ω < 1, =
1, and > 1 means negative purifying selection, neutral evolution, and positive selection.
"""
CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()
WRITE_CTL_FILE = 0
LOG_FILE = "paml_one_ratio_model.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def get_tree_path(trees_folder, species_folder_name):
    for tree in os.scandir(trees_folder):
        if tree.name.split('.')[0] == species_folder_name:
            return tree.name


def get_input_items(folder_in, trees_folder):
    """ parse root folder with files for paml
    parse tree_folder to get appropriate tree """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):
            tree_name = get_tree_path(trees_folder, species_folder.name)
            tree_path = os.path.join(trees_folder, tree_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.scandir(item):
                        if infile.name.split('.')[-1] == 'phy':  # and infile.name.split('.')[0].isnumeric():
                            yield folder_in, species_folder.name, item.name, infile.name, tree_path


def set_one_ratio_model(infile, phylo_tree, personal_dir):
    file_out_path = infile.replace('.phy', '_one_ratio.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylo_tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=1)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)  # 1:codons; 2:AAs; 3:codons-->AAs
    cml.set_options(CodonFreq=2)  # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    cml.set_options(clock=0)  # 0:no clock, 1:global clock; 2:local clock; 3:TipDate
    cml.set_options(aaDist=0)
    cml.set_options(model=0)  # models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches
    cml.set_options(NSsites=[0])  # * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
    # 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
    # 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
    # 13:3normal>0
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)  # 1: kappa fixed, 0: kappa to be estimated
    cml.set_options(kappa=2)  # initial or fixed kappa
    cml.set_options(fix_omega=0)  # 1:omega fixed, 0:omega to be estimated
    cml.set_options(omega=1)  # initial or fixed omega, for codons or codon-based AAs
    cml.set_options(fix_alpha=1)  # 0: estimate gamma shape parameter; 1: fix it at alpha
    cml.set_options(alpha=0.)
    cml.set_options(Malpha=0)
    cml.set_options(ncatG=8)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=1)
    cml.set_options(Small_Diff=.5e-6)
    cml.set_options(cleandata=0)
    cml.set_options(method=0)
    cml.write_ctl_file()
    return file_out_path


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except:
            frames = inspect.trace()
            argvalues = inspect.getargvalues(frames[0][0])
            gene_name = argvalues.locals['args'][0][2]
            logging.exception("gene {}".format(gene_name))
            logging.exception("exception {} in file {}/{}/{}".format(sys.exc_info()[0],
                                                                     argvalues.locals['args'][0][1],
                                                                     gene_name,
                                                                     argvalues.locals['args'][0][3]))
    return wrapped_func


@trace_unhandled_exceptions
def run_codeml(input_tuple, time_out, exec_path, overwrite_flag):
    folder_in, species_folder, item_folder, infile_name, phylogeny_tree_path = input_tuple
    item_folder_path = os.path.join(folder_in, species_folder, item_folder)
    file_gene_name = infile_name.split('.')[0]
    os.chdir(item_folder_path)
    logging.info("Working with {}/{}".format(species_folder, file_gene_name))
    infile_path = os.path.join(item_folder_path, infile_name)
    file_out_path = set_one_ratio_model(infile_path, phylogeny_tree_path, item_folder_path)

    if not overwrite_flag:
        if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
            with open(file_out_path, 'r') as o_f:
                if "Time used" in o_f.readlines()[-1]:
                    logging.info("The work has already been done for file {}...continue...".format(file_out_path))
                    return

    p = subprocess.Popen(exec_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    p.wait(timeout=time_out)  # Timeout in seconds
    # stdout, stderr = p.communicate()
    if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
        with open(file_out_path, 'r') as o_f:
            if "Time used" in o_f.readlines()[-1]:
                logging.info("OK gene {}".format(file_gene_name))
                logging.info("OK file {}/{}".format(species_folder, file_gene_name))
            else:
                logging.warning("FAIL for file {}/{}".format(species_folder, file_gene_name))
                raise codeml.CodemlError


def get_errors_from_log_file():
    print("get_errors_from_log_file")
    global ERROR_FILES_TO_WRITE
    with open(LOG_FILE, 'r') as lf:
        for line in lf:
            if 'ERROR : gene' in line:
                before_keyword, keyword, after_keyword = line.partition('ERROR : gene ')
                gene_name = after_keyword.replace('\n', '')
                ERROR_FILES_TO_WRITE.add(gene_name)
                print("add error", gene_name)


def get_corrects_from_log_file():
    print("get_corrects_from_log_file")
    global CORRECT_FILES_TO_WRITE
    with open(LOG_FILE, 'r') as lf:
        for line in lf:
            print('line', line)
            if 'INFO : OK gene' in line:
                before_keyword, keyword, after_keyword = line.partition('INFO : OK gene ')
                gene_name = after_keyword.replace('\n', '')
                CORRECT_FILES_TO_WRITE.add(gene_name)
                print("add correct", gene_name)


def write_correct_error_files(output_dir):
    get_errors_from_log_file()
    get_corrects_from_log_file()
    global CORRECT_FILES_TO_WRITE
    global ERROR_FILES_TO_WRITE
    logging.info("CORRECT_FILES_TO_WRITE {}, {}".format(len(CORRECT_FILES_TO_WRITE), CORRECT_FILES_TO_WRITE))
    logging.info("ERROR_FILES_TO_WRITE {}, {}".format(len(ERROR_FILES_TO_WRITE), ERROR_FILES_TO_WRITE))
    resulting_file = os.path.join(output_dir, 'paml_one_ratio_model_summary.xlsx')
    writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
    df_corr = pd.DataFrame({'Gene name': list(CORRECT_FILES_TO_WRITE)})
    df_corr.to_excel(writer, sheet_name='correct files', index=False)
    df_err = pd.DataFrame({'Gene name': list(ERROR_FILES_TO_WRITE)})
    df_err.to_excel(writer, sheet_name='exception files', index=False)
    logging.info("Summary has been written to {}".format(resulting_file))
    writer.save()


def main(folder_in, exec_path, trees_folder, time_out, threads_number, overwrite_flag):
    inputs = list(get_input_items(folder_in, trees_folder))  # list of tuples
    len_inputs = len(inputs)
    logging.info("Number of files should be analyzed = {}".format(len_inputs))
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    pool = multiprocessing.Pool(processes=threads_number)
    res = pool.starmap_async(run_codeml, zip(inputs, len_inputs * [time_out], len_inputs * [exec_path], len_inputs * [
        overwrite_flag]))
    res.wait()
    res.get()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', help='Path to the codeml executable', nargs='?', default='codeml')
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?', required=True)
    parser.add_argument('--tree', help='Path to the folder with trees for paml', nargs='?', required=True)
    parser.add_argument('--t', help='Timeout in seconds to get paml finish', nargs='?', default=120)
    parser.add_argument('--threads', help='Number of threads to use', nargs='?', required=True)
    parser.add_argument('--rework', help='"y" if overwrite existing files, default "n"', nargs='?', default='n')
    args = parser.parse_args()
    in_folder = args.i
    executable_path = args.e
    tree_folder = args.tree
    timeout = int(args.t)
    threads = int(args.threads)
    if args.rework == 'y':
        rework = True
    else:
        rework = False
    logging.info("Options in use:\n--e {} --i {} --tree {} --t {} --threads {} --rework {}".
                 format(executable_path, in_folder, tree_folder, timeout, threads, rework))
    try:
        main(in_folder, executable_path, tree_folder, timeout, threads, rework)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    write_correct_error_files(in_folder)
    logging.info("The work has been completed")
