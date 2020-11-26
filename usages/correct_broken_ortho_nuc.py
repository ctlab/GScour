#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

PROCESSED_FILES = dict()
BROKEN_SPECIES = list()
BROKEN_ACCORDANCE = dict()  # broken accordance with protein length (nuc = protein length * 3)
BROKEN_MULTIPLE_THREE = dict()
BROKEN_STOP_CODON = dict()
BROKEN_START_CODON = dict()
PREVIOUS_BROKEN = list()
NUMBER_OF_NEED_TO_BE_WRITTEN = 0
LOG_FILE = "get_ortho_nuc_from_broken.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)


def check_start_codon(seq, file_out_number, protein_id, start_codons="ATG"):
    global BROKEN_START_CODON
    codon = str(seq[0:3])
    if codon in start_codons:
        logging.info("ok start codon for file {} protein {}".format(file_out_number, protein_id))
        try:
            BROKEN_START_CODON.get(file_out_number).remove(protein_id)
            if not BROKEN_START_CODON.get(file_out_number):
                del BROKEN_START_CODON[file_out_number]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    if file_out_number not in BROKEN_START_CODON.keys():
        BROKEN_START_CODON[file_out_number] = list()
    logging.info("broken start codon for file {} protein {}".format(file_out_number, protein_id))
    if protein_id not in BROKEN_START_CODON.get(file_out_number):
        BROKEN_START_CODON.get(file_out_number).append(protein_id)


def check_accordance_with_protein_length(seq_length, protein_length, protein_id, file_out_number):
    global BROKEN_ACCORDANCE
    n = seq_length - 3 * protein_length
    if n == 0:
        logging.info(
            "check_accordance_with_protein_length-OK for file {} protein_id {}".format(file_out_number, protein_id))
        try:
            BROKEN_ACCORDANCE.get(file_out_number).remove(protein_id)
            if not BROKEN_ACCORDANCE.get(file_out_number):
                del BROKEN_ACCORDANCE[file_out_number]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    logging.warning("length of nucleotide_seq {} not equal proteins length*3={}*3 "
                    "with protein_id {} for file {}".
                    format(seq_length, protein_length, protein_id, file_out_number))
    if file_out_number not in BROKEN_ACCORDANCE.keys():
        BROKEN_ACCORDANCE[file_out_number] = list()
    if protein_id not in BROKEN_ACCORDANCE.get(file_out_number):
        BROKEN_ACCORDANCE.get(file_out_number).append(protein_id)


def check_stop_codon(seq, file_out_number, protein_id, stop_codons=("TAG", "TGA", "TAA")):
    global BROKEN_STOP_CODON
    codon = seq[-1 - 2:]
    if codon in stop_codons:
        logging.info("stop codon-OK for file {} protein_id {}".format(file_out_number, protein_id))
        try:
            BROKEN_STOP_CODON.get(file_out_number).remove(protein_id)
            if not BROKEN_STOP_CODON.get(file_out_number):
                del BROKEN_STOP_CODON[file_out_number]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    if file_out_number not in BROKEN_STOP_CODON.keys():
        BROKEN_STOP_CODON[file_out_number] = list()
    logging.info("broken stop codon for file {} protein {}".format(file_out_number, protein_id))
    if protein_id not in BROKEN_STOP_CODON.get(file_out_number):
        BROKEN_STOP_CODON.get(file_out_number).append(protein_id)


def check_multiple_of_three(seq, file_out_number, protein_id):
    global BROKEN_MULTIPLE_THREE
    if len(seq) % 3 == 0:
        logging.info("nucleotide_seq is multiple of three for file number {}".format(file_out_number))
        try:
            BROKEN_MULTIPLE_THREE.get(file_out_number).remove(protein_id)
            if not BROKEN_MULTIPLE_THREE.get(file_out_number):
                del BROKEN_MULTIPLE_THREE[file_out_number]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    if file_out_number not in BROKEN_MULTIPLE_THREE.keys():
        BROKEN_MULTIPLE_THREE[file_out_number] = list()
    if protein_id not in BROKEN_MULTIPLE_THREE.get(file_out_number):
        BROKEN_MULTIPLE_THREE.get(file_out_number).append(protein_id)
    logging.warning("length {} not multiple of three in file number {}:\n{}".
                    format(len(seq), file_out_number, seq))


def check_common_accordance(check_multiple, check_start, check_stop, check_accordance):
    if (check_start or not check_start) and (check_accordance or not check_accordance) and (
            check_stop or not check_stop) and check_multiple:
        logging.info("Check common accordance - OK: check_multiple - {}, check_start - {}, check_stop - {}, "
                     "check_accordance "
                     "- {}".format(check_multiple,
                                   check_start, check_stop, check_accordance))
        return True


def delete_from_broken(file_out_number, protein_id):
    global BROKEN_START_CODON
    global BROKEN_STOP_CODON
    global BROKEN_ACCORDANCE
    broken_list = [BROKEN_START_CODON, BROKEN_STOP_CODON, BROKEN_ACCORDANCE]
    try:
        for item in broken_list:
            item.get(file_out_number).remove(protein_id)
            if not item.get(file_out_number):
                del item[file_out_number]
    except (AttributeError, KeyError, ValueError):
        pass


def extract_seq_from_fasta(initfna_filepath, feature, record_id, species_numerating,
                           check_start, check_accordance, check_stop):
    fna_file_path = os.path.join(initfna_filepath, '{}.{}'.format(species_numerating, 'fna'))
    extracted_seq = ""
    for record in SeqIO.parse(fna_file_path, "fasta"):
        if record.id == record_id:
            for num, i in enumerate(feature.location.parts, start=0):
                if check_start and not check_accordance and check_stop:
                    logging.info("Manually extracted: check_start and not check_accordance and check_stop: check if "
                                 "in manually extracted "
                                 "sequence 'AATG' in the start and it is multiple of three")
                    extracted_seq += record.seq[i.start.position - 1:i.end.position]
                elif not check_start and not check_accordance and not check_stop:
                    extracted_seq += record.seq[i.start.position + 6:i.end.position + 1]
                    logging.info("Manually extracted: if not check_start and not check_accordance and not check_stop:"
                                 " all start pos + 6, all stop pos + 1")
                elif check_start and not check_accordance and not check_stop:
                    extracted_seq += record.seq[i.start.position:i.end.position + 1]
                    logging.info("Manually extracted: if check_start and not check_accordance and not check_stop:all "
                                 "start + 1, all stop + 1")
                else:
                    extracted_seq += record.seq[i.start.position - 1:i.end.position]
            logging.info("Manually extracted seq length = {}, seq is:\n{}".format(len(extracted_seq), extracted_seq))
            return extracted_seq


def check_translate(seq, protein_translation, initfna_filepath, feature, record_id,
                    species_numerating, check_multiple, check_start, check_stop, check_accordance, protein_length_gbff):
    logging.info("start checking translate with : check_multiple - {}, check_start - {},"
                 " check_stop - {}, check_accordance - {}".format(check_multiple, check_start, check_stop,
                                                                  check_accordance))
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
        'NNN': 'X',
        }
    protein = ""
    try:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table.get(codon)
            # if protein[j] != protein_translation[j]:
            # candidates = {k: v for k, v in table.items() if v == protein_translation[j]}
            # for k, v in candidates.items():
    except TypeError:
        pass

    protein_counted_length = len(protein)
    if protein_translation != protein:
        extracted_seq = extract_seq_from_fasta(initfna_filepath, feature, record_id, species_numerating,
                                               check_start, check_accordance, check_stop)
        extracted_seq_length = len(extracted_seq)
    else:
        logging.info("protein counted = protein translation from .gbff, return initial seq")
        return seq, True

    if check_stop:
        logging.info("from check_translate:\n protein counted length={}\nprotein translation length from .gbff={}\n "
                     "length seq[:-3]=length"
                     "without right stop codon(check_stop=True)={}\n"
                     "length of man extracted seq={}"
                     "protein counted:\n{}\n"
                     "protein translation:\n{}\n"
                     "nuc sequence, show stop codon:\n{}\n"
                     "manually extracted sequence:\n{}\n".
                     format(protein_counted_length, protein_length_gbff, len(seq[:-3]), extracted_seq_length, protein,
                            protein_translation,
                            seq, extracted_seq))
    else:
        logging.info("from check_translate:\n protein counted length={}\nprotein translation length from .gbff={}\n "
                     "length seq=length with broken stop codon(check_stop=False)={}\n"
                     "length of man extracted seq={}"
                     "protein counted:\n{}\n"
                     "protein translation:\n{}\n"
                     "nuc sequence, show broken stop codon\n{}"
                     "manually extracted sequence:\n{}\n".
                     format(protein_counted_length, protein_length_gbff, len(seq[:-3]), extracted_seq_length, protein,
                            protein_translation,
                            seq, extracted_seq))
    return extracted_seq, False


def anti_repeat_check(anti_repeat_store, protein_id, nucleotide_seq, file_out_number, start,
                      accordance, stop, multiple, gene, nucleotide_seq_length, file_number, species_numerating):
    if not anti_repeat_store.get(protein_id):
        anti_repeat_store[protein_id] = {
            "seq": nucleotide_seq, "start": start, "accordance": accordance,
            "stop": stop, "multiple": multiple, "gene": gene, "length": nucleotide_seq_length,
            "file_number": file_number, "species_numerating": species_numerating
            }
        return nucleotide_seq
    else:
        previous_seq = anti_repeat_store.get(protein_id).get("seq")
        if previous_seq == nucleotide_seq:
            logging.warning("Repeat file {} protein_id {} nucleotide_seqs are equal".format(
                file_out_number, protein_id))
            return nucleotide_seq
        else:
            logging.warning("Repeat file {} protein_id {} nucleotide_seqs are NOT equal:\n"
                            "length of previous = {}, length of current = {}\n"
                            "previous seq:\n{}\ncurrent seq:\n{}".format(file_out_number, protein_id, len(previous_seq),
                                                                         len(nucleotide_seq),
                                                                         previous_seq, nucleotide_seq))

            previous_start = anti_repeat_store.get(start)
            previous_accordance = anti_repeat_store.get(accordance)
            previous_stop = anti_repeat_store.get(stop)
            previous_multiple = anti_repeat_store.get(multiple)
            if [previous_start, previous_accordance, previous_stop, previous_multiple].count(True) > [start, accordance,
                                                                                                      stop, multiple]. \
                    count(True):
                logging.info("nucleotide seq was replaces by previous seq")
                return previous_seq
            else:
                anti_repeat_store[protein_id] = {
                    "seq": nucleotide_seq, "start": start, "accordance": accordance,
                    "stop": stop, "multiple": multiple, "gene": gene, "length": nucleotide_seq_length,
                    "file_number": file_number, "species_numerating": species_numerating
                    }
                return nucleotide_seq


def write_fasta_file(directory_out, file_out_number, seq):
    global PROCESSED_FILES
    with open(os.path.join(directory_out, file_out_number + ".fna"), "a") as ofile:
        SeqIO.write(seq, ofile, "fasta")
        if not WRITTEN_FILES.get(file_out_number):
            WRITTEN_FILES[file_out_number] = 0
        WRITTEN_FILES[file_out_number] += 1


def write_log_file(log_file, gene, protein_id, seq_length,
                   file_out_number, species_numerating):
    with open(log_file, "a") as ofile:
        file_size = os.path.getsize(log_file)
        if file_size == 0:
            ofile.write("gene_name - protein_id - seq_length - file number - species_numerating\n")
        ofile.write("{} - {} - {} - {} - {}\n".format(gene, protein_id, seq_length,
                                                      file_out_number, species_numerating))
        logging.info("wrote seq of gene {} corresponding {} in file number {}".format(gene, protein_id,
                                                                                      file_out_number))


def write_files(seq_store, directory_out):
    for protein_id in seq_store.keys():
        nucleotide_seq = seq_store.get(protein_id).get("seq")
        if seq_store.get(protein_id).get("stop"):
            nucleotide_seq = nucleotide_seq[:-3]
        species_numerating = seq_store.get(protein_id).get("species_numerating")
        gene = seq_store.get(protein_id).get("gene")
        file_number = seq_store.get(protein_id).get("file_number")
        nucleotide_seq_length = seq_store.get(protein_id).get("length")
        seq_record = SeqRecord(nucleotide_seq, id=species_numerating, description="")

        write_fasta_file(directory_out, file_number, seq_record)
        log_file = os.path.join(directory_out, file_number + ".log")
        write_log_file(log_file, gene, protein_id, nucleotide_seq_length,
                       file_number, species_numerating)


def get_and_write_nucleotide_seq(gb_file, ortho_protein_ids, directory_out, species_numerating,
                                 initfna_filepath):
    anti_repeat_store = dict()
    global BROKEN_START_CODON
    global BROKEN_STOP_CODON
    global BROKEN_MULTIPLE_THREE
    global BROKEN_ACCORDANCE
    global PREVIOUS_BROKEN
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("protein_id") in ortho_protein_ids:
                    protein_id = feature.qualifiers.get("protein_id")[0]  # str
                    protein_translation = feature.qualifiers.get("translation")[0]  # str
                    protein_translation_length = len(protein_translation)  # int
                    nucleotide_seq = feature.location.extract(record.seq)
                    index_count = np.where(ortho_protein_ids == protein_id)[0][0]
                    if int(index_count - 1) not in PREVIOUS_BROKEN:
                        continue
                    file_out_number = str(index_count + 1)
                    nucleotide_seq_length = len(nucleotide_seq)
                    gene = feature.qualifiers.get('gene')[0]

                    check_start = check_start_codon(nucleotide_seq, file_out_number, protein_id)
                    check_accordance = check_accordance_with_protein_length(nucleotide_seq_length - 3,
                                                                            protein_translation_length,
                                                                            protein_id, file_out_number)
                    check_stop = check_stop_codon(nucleotide_seq, file_out_number, protein_id)
                    check_multiple = check_multiple_of_three(nucleotide_seq, file_out_number, protein_id)
                    if check_common_accordance(check_multiple, check_start, check_stop, check_accordance):
                        delete_from_broken(file_out_number, protein_id)
                    else:
                        extracted_seq, check_trans_result = check_translate(nucleotide_seq, protein_translation,
                                                                            initfna_filepath, feature, record.id,
                                                                            species_numerating,
                                                                            check_multiple, check_start, check_stop,
                                                                            check_accordance,
                                                                            protein_translation_length)
                        if not check_trans_result:
                            logging.info("starting checks for manually extracted sequence")
                            check_start = check_start_codon(extracted_seq, file_out_number, protein_id)
                            check_accordance = check_accordance_with_protein_length(len(extracted_seq) - 3,
                                                                                    protein_translation_length,
                                                                                    protein_id, file_out_number)
                            check_stop = check_stop_codon(extracted_seq, file_out_number, protein_id)
                            check_multiple = check_multiple_of_three(extracted_seq, file_out_number, protein_id)
                            if check_common_accordance(check_multiple, check_start, check_stop, check_accordance):
                                delete_from_broken(file_out_number, protein_id)
                                nucleotide_seq = extracted_seq
                        """extracted_seq - No, check_multiple_of_three(nucleotide_seq) = False, -
                           nucleotide_seq - for further processing"""
                    # check duplicates
                    nucleotide_seq = anti_repeat_check(anti_repeat_store, protein_id,
                                                       nucleotide_seq, file_out_number,
                                                       check_start, check_accordance, check_stop, check_multiple, gene,
                                                       nucleotide_seq_length, file_out_number, species_numerating)

                    logging.info("Extracted and analysis has ended: detected gene {} corresponding protein_id {} with "
                                 "protein_length={}\n"
                                 "length of nucleotide seq={}\n"
                                 "seq:\n{}".format(gene, protein_id, protein_translation_length, nucleotide_seq_length,
                                                   nucleotide_seq))
    write_files(anti_repeat_store, directory_out)


def conv_int(val):
    if not val:
        return 0
    if val == 'Length':
        return 'Length'
    try:
        return np.int(val)
    except:
        return np.int(0)


def conv_string(val):
    if not val:
        return ""
    if val == 'Protein product':
        return 'Protein product'
    try:
        return np.str(val)
    except:
        return np.str("")


def replace_broken_files(directory_out):
    broken_species_folder = "new_broken_species_files"
    broken_multiple_folder = "new_broken_multiple_files"
    os.makedirs(broken_species_folder)
    os.makedirs(broken_multiple_folder)
    for file_number in BROKEN_SPECIES:
        os.replace(os.path.join(directory_out, file_number + ".fna"),
                   os.path.join(broken_species_folder, file_number + ".fna"))
        os.replace(os.path.join(directory_out, file_number + ".log"),
                   os.path.join(broken_species_folder, file_number + ".log"))
    for file_number in BROKEN_MULTIPLE_THREE:
        os.replace(os.path.join(directory_out, file_number + ".fna"),
                   os.path.join(broken_multiple_folder, file_number + ".fna"))
        os.replace(os.path.join(directory_out, file_number + ".log"),
                   os.path.join(broken_multiple_folder, file_number + ".log"))


def parse_dir(infolder):
    global PREVIOUS_BROKEN
    for infile in os.listdir(infolder):
        if infile.split('.')[-1] == 'fna':
            PREVIOUS_BROKEN.append(int(infile.split('.')[0]))


def main(orthodata_filepath, annotation_gbff, annotation_csv, initfna_filepath, species, directory_out):
    global NUMBER_OF_NEED_TO_BE_WRITTEN
    global PREVIOUS_BROKEN
    if not os.path.isdir(directory_out):
        os.makedirs(directory_out)
    ortho_data = pd.read_csv(orthodata_filepath, sep='\t', usecols=range(3, 3 + species))
    for column_number, _ in enumerate(ortho_data.columns):
        species_numerating = str(column_number + 1)
        annotation_gbff_path = os.path.join(annotation_gbff, "{}.{}".format(species_numerating, 'gbff'))
        annotation_csv_path = os.path.join(annotation_csv, "{}.{}".format(species_numerating, 'csv'))
        logging.info('Working with annotation files {}\n{}'.format(annotation_gbff_path, annotation_csv_path))
        logging.info('Working with species number {}'.format(species_numerating))
        df = pd.read_csv(annotation_csv_path, sep=',', names=['Protein product', 'Length'], usecols=[8, 9],
                         converters={'Protein product': conv_string, 'Length': conv_int}, low_memory=False)
        ortho_protein_ids = ortho_data.iloc[:, column_number].values

        NUMBER_OF_NEED_TO_BE_WRITTEN = len(PREVIOUS_BROKEN)
        get_and_write_nucleotide_seq(annotation_gbff_path, ortho_protein_ids, directory_out,
                                     species_numerating, initfna_filepath)

    logging.info("NUMBER_OF_NEED_TO_BE_WRITTEN: {}".format(NUMBER_OF_NEED_TO_BE_WRITTEN))

    for file, written_species in PROCESSED_FILES.items():
        if written_species != species:
            if file not in BROKEN_SPECIES:
                BROKEN_SPECIES.append(file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Path to the single_copy_orthologs.tsv', nargs='?')
    parser.add_argument('--gbff', help='Path to the folder with annotation .gbff files from '
                                       'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--csv', help='Path to the folder with annotation .csv files from '
                                      'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--genome', help='Path to the folder with reference genome'
                                         'in FASTA format', nargs='?')
    parser.add_argument('--species', help='Number of species', nargs='?')
    parser.add_argument('--out', help='Path to the folder for result write out', nargs='?')
    parser.add_argument('--broken', help='Path to the folder with broken files of previous launch', nargs='?',
                        default=PREVIOUS_BROKEN)
    args = parser.parse_args()

    try:
        if type(args.broken) != list:
            parse_dir(args.broken)
        logging.info("previous_broken list of length {}:\n{}".format(len(PREVIOUS_BROKEN), PREVIOUS_BROKEN))
        main(args.ortho, args.gbff, args.csv, args.genome, int(args.species), args.out)
        written_files_number = len(PROCESSED_FILES)
        delta = NUMBER_OF_NEED_TO_BE_WRITTEN - written_files_number
        if delta == 0:
            logging.info("All files are written")

        logging.info("NUMBER_OF_NEED_TO_BE_WRITTEN = {},  WRITTEN_FILES = {}, where {} in BROKEN_SPECIES list: {}"
                     .format(NUMBER_OF_NEED_TO_BE_WRITTEN, written_files_number, len(BROKEN_SPECIES),
                             repr(BROKEN_SPECIES)))
        if BROKEN_SPECIES or BROKEN_MULTIPLE_THREE:
            replace_broken_files(args.out)
            residue = written_files_number - len(BROKEN_SPECIES) - len(BROKEN_MULTIPLE_THREE)
            logging.info("removed broken species, broken multiple files into folders 'broken_species_files'"
                         "'broken_multiple_files' in cwd,"
                         "please check out folder for .fna files number: {}".format(residue))
    except:
        logging.exception("Unexpected error")

    logging.warning("BROKEN_SPECIES {} : {}".format(len(BROKEN_SPECIES), BROKEN_SPECIES))
    logging.warning("BROKEN_STOP_CODON {} : {}".format(len(BROKEN_STOP_CODON), BROKEN_STOP_CODON))
    logging.warning("BROKEN_START_CODON {} : {}".format(len(BROKEN_START_CODON), BROKEN_START_CODON))
    logging.warning("BROKEN_ACCORDANCE {} : {}".format(len(BROKEN_ACCORDANCE), BROKEN_ACCORDANCE))
    logging.warning("BROKEN_MULTIPLE_THREE {} : {}".format(len(BROKEN_MULTIPLE_THREE), BROKEN_MULTIPLE_THREE))
    logging.info("WRITTEN_FILES {}:\n{}".format(len(PROCESSED_FILES), repr(PROCESSED_FILES)))
    logging.info("The work has been completed")
