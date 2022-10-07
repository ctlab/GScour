#!/usr/bin/env python
import argparse
import subprocess
import functools, multiprocessing, inspect
import sys
import pandas as pd
from Bio.Phylo.PAML import codeml
import os
import logging

CORRECT_FILES_TO_WRITE = set()
ERROR_FILES_TO_WRITE = set()

LOG_FILE = "paml_branch_site_model.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)

"""There are two hypothesis:
H0 (The null model for Branch-site model A): 
    Model A1: model = 2, NSsites = 2, fix_omega = 1, omega = 1
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
H1 (Alternative model, Model A: model = 2, NSsites = 2, fix_omega = 0 ): 
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
    
From readme of paml example lysozymeLarge.ctl:
Alternative hypothesis (branch site model A, with w2 estimated):
model = 2    NSsites = 2   fix_omega = 0   omega = 1.5 (or any value > 1)

As the branch-site model is
known to cause computational difficulties for the numer-
ical optimization algorithm in PAML, each analysis is con-
ducted three times with different starting values to ensure
that the global peak is found (Statistical Properties of the Branch-Site Test of Positive
Selection, Ziheng Yang and Mario dos Reis)
"""


def get_tree_path(trees_folder, species_folder_name):
    for tree in os.scandir(trees_folder):
        if tree.name.split('.')[0] == species_folder_name:
            return tree.name


def get_input_items(folder_in, trees_folder):
    """ parse root folder with files for paml
    parse tree_folder to get appropriate tree """
    for species_folder in os.scandir(folder_in):
        if os.path.isdir(species_folder):  # and species_folder.name.isdigit():
            logging.info("working with species folder {}".format(species_folder.name))
            tree_name = get_tree_path(trees_folder, species_folder.name)
            tree_path = os.path.join(trees_folder, tree_name)
            for item in os.scandir(species_folder):
                if os.path.isdir(item):
                    for infile in os.listdir(item):
                        if infile.split('.')[-1] == 'phy':
                            yield folder_in, species_folder.name, item.name, infile, tree_path


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except:
            frames = inspect.trace()
            argvalues = inspect.getargvalues(frames[0][0])
            gene_name = argvalues.locals['args'][0][2]
            logging.exception("FAIL gene {}".format(gene_name))
            logging.exception("exception {} in file {}/{}/{}".format(sys.exc_info()[0],
                                                                     argvalues.locals['args'][0][1],
                                                                     gene_name,
                                                                     argvalues.locals['args'][0][3]))

    return wrapped_func


def set_alternative_hypothesis(infile, phylo_tree, personal_dir):
    file_out_path = infile.replace('.phy', '_alter1.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylo_tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=0)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=2)
    cml.set_options(NSsites=[2])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=0)
    cml.set_options(omega=16)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    cml.write_ctl_file()
    return cml, file_out_path


def set_null_hypothesis(infile, phylo_tree, personal_dir):
    file_out_path = infile.replace('.phy', '_null1.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=phylo_tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=0)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=2)
    cml.set_options(NSsites=[2])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=1)
    cml.set_options(omega=1)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    cml.write_ctl_file()
    return cml, file_out_path


@trace_unhandled_exceptions
def run_paml(input_tuple, exec_path, hypothesis_type, overwrite_flag, time_out):
    folder_in, species_folder, item_folder, infile_name, phylogeny_tree_path = input_tuple
    item_folder_path = os.path.join(folder_in, species_folder, item_folder)
    file_gene_name = item_folder
    os.chdir(item_folder_path)
    file_id = "{}/{}".format(species_folder, file_gene_name)
    logging.info("Working with {}".format(file_id))
    infile_path = os.path.join(item_folder_path, infile_name)

    if hypothesis_type == "null":
        cml, file_out_path = set_null_hypothesis(infile_path, phylogeny_tree_path, item_folder_path)
    elif hypothesis_type == "alter":
        cml, file_out_path = set_alternative_hypothesis(infile_path, phylogeny_tree_path, item_folder_path)
    else:
        logging.warning("Check the type of hypothesis: null, alter")
        return

    if not overwrite_flag:
        if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
            with open(file_out_path, 'r') as o_f:
                if "Time used" in o_f.readlines()[-1]:
                    logging.info("The work has already been done for file {}...continue...".format(file_out_path))
                    return
    try:
        p = subprocess.Popen(exec_path, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        # return_code = p.wait(timeout=time_out)  # Timeout in seconds
        p.communicate(input=b'\n', timeout=time_out)
        if os.path.isfile(file_out_path) and os.path.getsize(file_out_path) > 0:
            with open(file_out_path, 'r') as o_f:
                if "Time used" in o_f.readlines()[-1]:
                    logging.info("OK gene {}".format(file_gene_name))
                    logging.info("OK paml for file {}".
                                 format(file_id))

                else:
                    logging.error("FAIL gene {}".format(file_gene_name))
                    logging.error("FAIL paml for file {}".format(file_id))

    except subprocess.SubprocessError:
        # logging.warning("FAIL gene {}".format(file_gene_name))
        logging.warning("SubprocessError for {} hypothesis_type {}, RESUMING...".format(file_id, hypothesis_type))
        # ":\n{
        # }\nerr.args:\n{
        # }\nerr.stderr:\n{}".
        # format(file_id, hypothesis_type, err, err.stderr, err.args))
    # except subprocess.CalledProcessError as err:
    #     logging.exception("FAIL gene {}".format(file_gene_name))
    #     logging.exception("CalledProcessError for {} hypothesis_type {}:\n{}\nerr.args:\n{}\nerr.stderr:\n{}".
    #                       format(file_id, hypothesis_type, err, err.stderr, err.args))
    # except subprocess.TimeoutExpired as err:
    #     p.kill()
    #     logging.exception("FAIL gene {}".format(file_gene_name))
    #     logging.exception("Killed {} hypothesis_type {}, err.args: {}".
    #                       format(file_id, hypothesis_type, err.args))


def get_correct_from_log_file():
    global CORRECT_FILES_TO_WRITE
    with open(LOG_FILE, 'r') as lf:
        for line in lf:
            if 'OK gene ' in line:
                before_keyword, keyword, after_keyword = line.partition('OK gene ')
                gene_name = after_keyword.replace('\n', '')
                CORRECT_FILES_TO_WRITE.add(gene_name)


def get_errors_from_log_file():
    global ERROR_FILES_TO_WRITE
    global CORRECT_FILES_TO_WRITE
    with open(LOG_FILE, 'r') as lf:
        for line in lf:
            if 'FAIL gene ' in line:
                before_keyword, keyword, after_keyword = line.partition('FAIL gene ')
                gene_name = after_keyword.replace('\n', '')
                ERROR_FILES_TO_WRITE.add(gene_name)
                if gene_name in CORRECT_FILES_TO_WRITE:
                    CORRECT_FILES_TO_WRITE.remove(gene_name)


def write_correct_and_error_files(output_dir):
    get_correct_from_log_file()
    get_errors_from_log_file()
    global CORRECT_FILES_TO_WRITE
    global ERROR_FILES_TO_WRITE
    logging.info("CORRECT_FILES_TO_WRITE {}".format(len(CORRECT_FILES_TO_WRITE)))
    logging.info("ERROR_FILES_TO_WRITE {}".format(len(ERROR_FILES_TO_WRITE)))
    resulting_file = os.path.join(output_dir, 'paml_branch_site_model_summary.xlsx')
    writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
    df_corr = pd.DataFrame({'Gene name': list(CORRECT_FILES_TO_WRITE)})
    df_corr.to_excel(writer, sheet_name='correct files', index=False)
    df_err = pd.DataFrame({'Gene name': list(ERROR_FILES_TO_WRITE)})
    df_err.to_excel(writer, sheet_name='exception files', index=False)
    logging.info("Summary has been written to {}".format(resulting_file))
    writer.save()


def main(folder_in, trees_folder, exec_path, number_of_threads, overwrite_flag, time_out):
    inputs = list(get_input_items(folder_in, trees_folder))
    len_inputs = len(inputs)
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)

    pool = multiprocessing.Pool(processes=number_of_threads)
    i = pool.starmap_async(run_paml, zip(inputs, len_inputs * [exec_path],
                                         len_inputs * ["null"], len_inputs * [overwrite_flag], len_inputs * [time_out]))
    i.wait()
    i.get()
    i = pool.starmap_async(run_paml, zip(inputs, len_inputs * [exec_path],
                                         len_inputs * ["alter"], len_inputs * [overwrite_flag],
                                         len_inputs * [time_out]))
    i.wait()
    i.get()
    logging.info("Number of files for analysis: {}".format(len_inputs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', help='Path to the codeml executable', nargs='?', default="codeml")
    parser.add_argument('--timeout', help='Timeout for codeml in seconds, default=500', nargs='?', default='500')
    parser.add_argument('--i', help='The full path to the folder contains folders with input files for paml',
                        nargs='?', required=True)
    parser.add_argument('--tree', help='Path to the folder with trees for paml', nargs='?', required=True)
    parser.add_argument('--threads', help='Number of threads to use', nargs='?', required=True)
    parser.add_argument('--rework', help='"y" if overwrite existing files, default "n"', nargs='?', default='n')
    args = parser.parse_args()
    in_folder = args.i
    executable_path = args.e
    tree_folder = args.tree
    threads = int(args.threads)
    timeout = int(args.timeout)
    if args.rework == 'y':
        rework = True
    else:
        rework = False
    logging.info("Path to the folder with input files for paml: {}\nPath to the tree: {}\nExecutable path: {}\n"
                 "Threads to use = {}, rework = {}, timeout={}".
                 format(in_folder, tree_folder, executable_path, threads, rework, timeout))
    try:
        main(in_folder, tree_folder, executable_path, threads, rework, timeout)
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
        raise e
    write_correct_and_error_files(in_folder)
    logging.info("The work has been completed")
