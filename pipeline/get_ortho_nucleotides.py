#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

CORRECT_FILES = dict()  # store the number of writings (1 writing = 1 seq = 1 species):
# for species=5 and group=5 there should be PROCESSED_FILES[file_name] = ['1,2,3,4,5', 5] (without BROKEN lists)
BROKEN_SPECIES = list()
CORRECT_RECORD = list()
BROKEN_ACCORDANCE = dict()  # broken accordance with protein length (nuc = protein length * 3)
BROKEN_MULTIPLE_THREE = dict()
BROKEN_STOP_CODON = dict()
BROKEN_START_CODON = dict()
ABSENT_IN_CDS = dict()
NUMBER_OF_NEED_TO_BE_WRITTEN = 0
LOG_FILE = "get_ortho_nucleotides.log"
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO, filename=LOG_FILE)
BROKEN_LIST = []  # broken_list for re-extracting seqs from list ['item1', 'item2'..], e.g ['3557', '5781', '1503']


def check_start_codon(seq, file_out_name, protein_id, start_codons="ATG"):
    global BROKEN_START_CODON
    codon = str(seq[0:3])
    if codon in start_codons:
        logging.info("ok start codon for file {} protein {}".format(file_out_name, protein_id))
        try:
            BROKEN_START_CODON.get(file_out_name).remove(protein_id)
            if not BROKEN_START_CODON.get(file_out_name):
                del BROKEN_START_CODON[file_out_name]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    if file_out_name not in BROKEN_START_CODON.keys():
        BROKEN_START_CODON[file_out_name] = list()
    logging.info("broken start codon for file {} protein {}".format(file_out_name, protein_id))
    if protein_id not in BROKEN_START_CODON.get(file_out_name):
        BROKEN_START_CODON.get(file_out_name).append(protein_id)


def check_accordance_with_protein_length(seq_length, protein_length, protein_id, file_out_name):
    global BROKEN_ACCORDANCE
    n = seq_length - 3 * protein_length
    if n in range(-9, 10):  # do softer condition to exclude very rough mistaken sequences
        if n == 0:
            logging.info(
                "check_accordance_with_protein_length-OK for file {} protein_id {}".format(file_out_name, protein_id))
        else:
            logging.info("check_accordance_with_protein_length-NEAR OK: delta = seq_length - 3 * protein_length = {}\n"
                         "for file {} protein_id {}".format(n, file_out_name, protein_id))
        try:
            BROKEN_ACCORDANCE.get(file_out_name).remove(protein_id)
            if not BROKEN_ACCORDANCE.get(file_out_name):
                del BROKEN_ACCORDANCE[file_out_name]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    logging.warning("length of nucleotide_seq {} not equal proteins length*3={}*3 "
                    "with protein_id {} for file {}".
                    format(seq_length, protein_length, protein_id, file_out_name))
    if file_out_name not in BROKEN_ACCORDANCE.keys():
        BROKEN_ACCORDANCE[file_out_name] = list()
    if protein_id not in BROKEN_ACCORDANCE.get(file_out_name):
        BROKEN_ACCORDANCE.get(file_out_name).append(protein_id)


def check_stop_codon(seq, file_out_name, protein_id, stop_codons=("TAG", "TGA", "TAA")):
    global BROKEN_STOP_CODON
    codon = seq[-1 - 2:]
    if codon in stop_codons:
        logging.info("stop codon-OK for file {} protein_id {}".format(file_out_name, protein_id))
        try:
            BROKEN_STOP_CODON.get(file_out_name).remove(protein_id)
            if not BROKEN_STOP_CODON.get(file_out_name):
                del BROKEN_STOP_CODON[file_out_name]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    if file_out_name not in BROKEN_STOP_CODON.keys():
        BROKEN_STOP_CODON[file_out_name] = list()
    logging.info("broken stop codon for file {} protein {}".format(file_out_name, protein_id))
    if protein_id not in BROKEN_STOP_CODON.get(file_out_name):
        BROKEN_STOP_CODON.get(file_out_name).append(protein_id)


def check_multiple_of_three(seq, file_out_name, protein_id):
    global BROKEN_MULTIPLE_THREE
    if len(seq) % 3 == 0:
        logging.info("nucleotide_seq is multiple of three for file {}".format(file_out_name))
        try:
            BROKEN_MULTIPLE_THREE.get(file_out_name).remove(protein_id)
            if not BROKEN_MULTIPLE_THREE.get(file_out_name):
                del BROKEN_MULTIPLE_THREE[file_out_name]
        except (AttributeError, KeyError, ValueError):
            pass
        return True
    if file_out_name not in BROKEN_MULTIPLE_THREE.keys():
        BROKEN_MULTIPLE_THREE[file_out_name] = list()
    if protein_id not in BROKEN_MULTIPLE_THREE.get(file_out_name):
        BROKEN_MULTIPLE_THREE.get(file_out_name).append(protein_id)
    logging.warning("length {} not multiple of three in file {}:\n{}".
                    format(len(seq), file_out_name, seq))


def check_common_accordance(check_multiple, check_start, check_stop, check_accordance):
    if (check_start or not check_start) and (check_stop or not check_stop) and check_multiple and check_accordance:
        logging.info("Check common accordance - OK: check_multiple - {}, check_start - {}, check_stop - {}, "
                     "check_accordance "
                     "- {}".format(check_multiple,
                                   check_start, check_stop, check_accordance))
        return True


def delete_from_broken(file_out_name, protein_id):
    global BROKEN_START_CODON
    global BROKEN_STOP_CODON
    global BROKEN_ACCORDANCE
    broken_list = [BROKEN_START_CODON, BROKEN_STOP_CODON, BROKEN_ACCORDANCE]
    try:
        for item in broken_list:
            item.get(file_out_name).remove(protein_id)
            if not item.get(file_out_name):
                del item[file_out_name]
    except (AttributeError, KeyError, ValueError):
        pass


def extract_seq_from_fasta(fna_file_path, feature, record_id, check_start, check_accordance, check_stop):
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


def check_translate(seq, protein_translation, fna_filepath, feature, record_id,
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
        fna_file_path = os.path.join(fna_filepath, '{}.{}'.format(species_numerating, 'fna'))
        if not os.path.isfile(fna_file_path):
            logging.info("From check_translate: File path {} does not exist\nReturn initial seq from .gbff".format(
                fna_file_path))
            return seq, False
        extracted_seq = extract_seq_from_fasta(fna_file_path, feature, record_id, check_start, check_accordance,
                                               check_stop)
        if not extracted_seq:
            logging.info("From check_translate: extract_seq_from_fasta failed\nReturn initial seq from .gbff".format(
                fna_file_path))
            return seq, False
        extracted_seq_length = len(extracted_seq)
    else:
        logging.info("protein counted = protein translation from .gbff, return initial seq")
        return seq, True

    if check_stop:
        logging.info("from check_translate:\n protein counted length={}\nprotein translation length from .gbff={}\n "
                     "length seq[:-3]=length"
                     "without right stop codon(check_stop=True)={}\n"
                     "length of man extracted seq={}\n"
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
                     "length of man extracted seq={}\n"
                     "protein counted:\n{}\n"
                     "protein translation:\n{}\n"
                     "nuc sequence, show broken stop codon\n{}"
                     "manually extracted sequence:\n{}\n".
                     format(protein_counted_length, protein_length_gbff, len(seq[:-3]), extracted_seq_length, protein,
                            protein_translation,
                            seq, extracted_seq))
    return extracted_seq, False


def anti_repeat_check(protein_id, seq_record, anti_repeat_store):
    if not anti_repeat_store.get(protein_id):
        anti_repeat_store[protein_id] = {
            "seq": seq_record  # seq_record = SeqRecord(record.seq, id=species_numerating, description="")
            }
        return True
    return False


def anti_repeat_check_with_comparing(seq_store, protein_id, nucleotide_seq, file_out_name, start,
                                     accordance, stop, multiple, gene, nucleotide_seq_length, species_numerating):
    if not seq_store.get(protein_id):
        seq_store[protein_id] = {
            "seq": nucleotide_seq, "start": start, "accordance": accordance,
            "stop": stop, "multiple": multiple, "gene": gene, "length": nucleotide_seq_length,
            "file_name": file_out_name, "species_numerating": species_numerating
            }
        return nucleotide_seq
    else:
        previous_seq = seq_store.get(protein_id).get("seq")
        if previous_seq == nucleotide_seq:
            logging.warning("Anti-repeat check: repeating file {} protein_id {} nucleotide_seqs are equal".format(
                file_out_name, protein_id))
            return nucleotide_seq
        else:
            logging.warning("Anti-repeat check: repeating file {} protein_id {} nucleotide_seqs are NOT equal:\n"
                            "length of previous = {}, length of current = {}\n"
                            "previous seq:\n{}\ncurrent seq:\n{}".format(file_out_name, protein_id, len(previous_seq),
                                                                         len(nucleotide_seq),
                                                                         previous_seq, nucleotide_seq))

            previous_start = seq_store.get(protein_id).get(start)
            previous_accordance = seq_store.get(protein_id).get(accordance)
            previous_stop = seq_store.get(protein_id).get(stop)
            previous_multiple = seq_store.get(protein_id).get(multiple)
            if [previous_start, previous_accordance, previous_stop, previous_multiple].count(True) > [start, accordance,
                                                                                                      stop, multiple]. \
                    count(True):
                logging.info("Anti-repeat check: replacing the current seq with the previous one: return previous "
                             "seq")
                return previous_seq
            else:
                seq_store[protein_id] = {
                    "seq": nucleotide_seq, "start": start, "accordance": accordance,
                    "stop": stop, "multiple": multiple, "gene": gene, "length": nucleotide_seq_length,
                    "file_name": file_out_name, "species_numerating": species_numerating
                    }
                logging.info("Anti-repeat check: current seq is better then previous, replacing the previous seq with "
                             "the current one: return current seq")
                return nucleotide_seq


def write_fasta_file(dir_out, file_out_name, seq_record):
    global CORRECT_FILES
    with open(os.path.join(dir_out, file_out_name + ".fna"), "a") as f:
        SeqIO.write(seq_record, f, "fasta")
        logging.info("Sequence corresponding to the species number {} has been added to the file {}".format(
            seq_record.id, file_out_name))
        if not CORRECT_FILES.get(file_out_name):
            CORRECT_FILES[file_out_name] = ["", 0]
        CORRECT_FILES[file_out_name][0] += '{},'.format(seq_record.id)
        CORRECT_FILES[file_out_name][1] += 1


def write_log_file(log_file, protein_id, seq_length, file_out_name, species_numerating):
    with open(log_file, "a") as o_file:
        file_size = os.path.getsize(log_file)
        if file_size == 0:
            o_file.write("gene_name - protein_id - seq_length - file_name - species_numerating\n")
        o_file.write("{} - {} - {} - {} - {}\n".format(file_out_name, protein_id, seq_length,
                                                       file_out_name, species_numerating))
        logging.info("Log file for gene-protein {}-{} recorded".format(file_out_name, protein_id))


def write_fasta_and_log(seq_store, protein_id, dir_out):
    nucleotide_seq = seq_store.get(protein_id).get("seq")
    if seq_store.get(protein_id).get("stop"):
        nucleotide_seq = nucleotide_seq[:-3]
    species_numerating = seq_store.get(protein_id).get("species_numerating")
    # gene = seq_store.get(protein_id).get("gene")
    file_name = seq_store.get(protein_id).get("file_name")
    nucleotide_seq_length = seq_store.get(protein_id).get("length")
    seq_record = SeqRecord(nucleotide_seq, id=species_numerating, description='')

    write_fasta_file(dir_out, file_name, seq_record)
    log_file = os.path.join(dir_out, file_name + ".log")
    write_log_file(log_file, protein_id, nucleotide_seq_length, file_name, species_numerating)


def get_seq_from_gbff(gb_file, ortho_protein_ids, ortho_gene_names):
    if not os.path.isfile(gb_file):
        logging.info("There is no such file path: {}\nReturn".format(gb_file))
        return
    logging.info("Starting to get sequences from .gbff {}".format(gb_file))
    df_ortho = pd.DataFrame({'species': ortho_protein_ids[:], 'Gene name': ortho_gene_names[:]})
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("protein_id") in df_ortho.values:
                    protein_id = feature.qualifiers.get("protein_id")[0]
                    gene_name = (df_ortho.loc[df_ortho['species'] == protein_id]['Gene name']).values[0]
                    protein_translation = feature.qualifiers.get("translation")[0]  # str
                    protein_translation_length = len(protein_translation)  # int
                    nucleotide_seq = feature.location.extract(record.seq)
                        # index_count = np.where(ortho_protein_ids == protein_id)[0][0]
                        # file_out_number = str(index_count + 1)
                    file_out_name = gene_name
                    nucleotide_seq_length = len(nucleotide_seq)
                        # gene = feature.qualifiers.get('gene')[0]
                    logging.info("get_seq_from_gbff: protein_id {} gene {}".format(protein_id, gene_name))
                    yield (protein_id, protein_translation, protein_translation_length, nucleotide_seq,
                           file_out_name, nucleotide_seq_length, gene_name, feature, record.id)


def get_seq_record_from_cds(cds_from_genomic_file, protein_id, species_numerating):
    seq_record, gene_name = "", ""
    for record in SeqIO.parse(cds_from_genomic_file, "fasta"):
        if re.search(r'_({})_'.format(protein_id), record.name):
            logging.info("get_seq_record_from_cds: protein_id {}, record {}".format(protein_id, record.name))
            # if protein_id in record.name:
            if "gene=" in record.description:
                gene_name = re.search(r'gene=([A-Za-z0-9]+)', record.description).group(1)
                logging.info("get_seq_record_from_cds: protein_id {} gene {}".format(protein_id, gene_name))
            else:
                logging.info("get_seq_record_from_cds protein_id {}, no gene name".format(protein_id))
            seq_record = SeqRecord(record.seq, id=species_numerating, description='')
            break
    return seq_record


def get_from_cds_and_write(cds_from_genomic_file, ortho_protein_ids, ortho_gene_names, species_numerating, dir_out,
                           seq_store):
    global ABSENT_IN_CDS
    if os.path.isfile(cds_from_genomic_file):
        logging.info("Starting to get sequences from cds_from_genomic {}".format(cds_from_genomic_file))
        for protein_id, gene_name in zip(ortho_protein_ids, ortho_gene_names):
            if protein_id == '*' or not protein_id or not isinstance(protein_id, str):
                continue
            seq_record = get_seq_record_from_cds(cds_from_genomic_file, protein_id, species_numerating)
            file_out_name = gene_name
            if not seq_record:  # no seq_record when this species is not in the group
                if not ABSENT_IN_CDS.get(species_numerating):
                    ABSENT_IN_CDS[species_numerating] = list()
                ABSENT_IN_CDS.get(species_numerating).append(protein_id)
                logging.info("protein_id {} is absent in {}".format(protein_id, cds_from_genomic_file))
                continue
            else:
                seq_length = len(seq_record.seq)
            # if not seq_record.description:  # gene_name
            #     logging.warning("empty gene_name from cds_from_genomic")
            if anti_repeat_check(protein_id, seq_record, seq_store):
                write_fasta_file(dir_out, file_out_name, seq_record)
                log_file = os.path.join(dir_out, file_out_name + ".log")
                write_log_file(log_file, protein_id, seq_length, file_out_name, species_numerating)
        logging.info("all sequences for species {} recorded".format(species_numerating))
        return True
    logging.info("No such cds_from_genomic {}".format(cds_from_genomic_file))
    return False


def get_and_write_nucleotide_seq(gb_file, cds_from_genomic_file, ortho_protein_ids, ortho_gene_names, dir_out,
                                 species_numerating,
                                 genome_fna_path):
    """
    :param ortho_gene_names:
    :param species_numerating:
    :param dir_out:
    :param ortho_protein_ids:
    :param gb_file: GenBank annotation file .gbff
    :param cds_from_genomic_file: RefSeq annotated cds_from_genomic.fna (or other source FASTA format of'
                                      'the nucleotide sequences corresponding to all CDS'
                                      'features annotated on the assembly)
    :param genome_fna_path: folder with genomes
    1. Trying to get nucleotide sequence from cds_from_genomic.fna
    2. If there is no such file, trying to find sequence by protein_id in .gbff
    3. Checking nucleotide sequence for start/stop codon/multiple 3/accordance with protein_length extracted from .gbff
    4. If checking failed, try to check with manually translating sequence
    5. If translating check failed, try to extract_seq_from_fasta (initial genome .fna file)
    """
    global BROKEN_START_CODON
    global BROKEN_STOP_CODON
    global BROKEN_MULTIPLE_THREE
    global BROKEN_ACCORDANCE
    seq_store = dict()
    """ trying to get nucleotide sequence from cds_from_genomic.fna """
    if get_from_cds_and_write(cds_from_genomic_file, ortho_protein_ids, ortho_gene_names, species_numerating,
                              dir_out, seq_store):
        return
    """ trying to find sequence by protein_id in .gbff """
    for (protein_id, protein_translation, protein_translation_length, nucleotide_seq,
         file_out_name, nucleotide_seq_length, gene, feature, record_id) in \
            get_seq_from_gbff(gb_file, ortho_protein_ids, ortho_gene_names):
        """ start checking """
        check_start = check_start_codon(nucleotide_seq, file_out_name, protein_id)
        check_accordance = check_accordance_with_protein_length(nucleotide_seq_length - 3,
                                                                protein_translation_length,
                                                                protein_id, file_out_name)
        check_stop = check_stop_codon(nucleotide_seq, file_out_name, protein_id)
        check_multiple = check_multiple_of_three(nucleotide_seq, file_out_name, protein_id)
        if check_common_accordance(check_multiple, check_start, check_stop, check_accordance):
            delete_from_broken(file_out_name, protein_id)  # common: if check_multiple, check_accordance-ok, then OK
        else:  # extract from genome .fna
            extracted_seq, check_trans_result = check_translate(nucleotide_seq, protein_translation,
                                                                genome_fna_path, feature, record_id,
                                                                species_numerating,
                                                                check_multiple, check_start, check_stop,
                                                                check_accordance,
                                                                protein_translation_length)
            if not check_trans_result and extracted_seq != nucleotide_seq:
                logging.info("starting checks for manually extracted sequence")
                check_start = check_start_codon(extracted_seq, file_out_name, protein_id)
                check_accordance = check_accordance_with_protein_length(len(extracted_seq) - 3,
                                                                        protein_translation_length,
                                                                        protein_id, file_out_name)
                check_stop = check_stop_codon(extracted_seq, file_out_name, protein_id)
                check_multiple = check_multiple_of_three(extracted_seq, file_out_name, protein_id)
                if check_common_accordance(check_multiple, check_start, check_stop, check_accordance):
                    delete_from_broken(file_out_name, protein_id)
                    nucleotide_seq = extracted_seq
                """ else: if the attempt to extract right seq from genome .fna has failed - extracted_seq - return to
                   nucleotide_seq from .gbff """
        # check duplicates
        nucleotide_seq = anti_repeat_check_with_comparing(seq_store, protein_id, nucleotide_seq, file_out_name,
                                                          check_start, check_accordance, check_stop, check_multiple,
                                                          gene, nucleotide_seq_length, species_numerating)

        logging.info("Extraction and analysis complete: detected gene {} corresponding protein_id {} with "
                     "protein_length {} "
                     "nucleotide seq of length {}\n"
                     "seq:\n{}\n".format(gene, protein_id, protein_translation_length,
                                         nucleotide_seq_length, nucleotide_seq))
        write_fasta_and_log(seq_store, protein_id, dir_out)
        seq_store = dict()


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


def replace_broken_files_and_write_table(dir_out, broken_folder, result_table_writer):
    for file_name in BROKEN_SPECIES:
        os.replace(os.path.join(dir_out, file_name + ".fna"),
                   os.path.join(broken_folder, file_name + ".fna"))
        os.replace(os.path.join(dir_out, file_name + ".log"),
                   os.path.join(broken_folder, file_name + ".log"))
    res = pd.DataFrame({'Gene name': BROKEN_SPECIES})
    res.to_excel(result_table_writer, sheet_name='broken species number', index=False)
    # res.to_csv(broken_species_file, sep="\t")


def main(orthodata_filepath, annotation_gbff, cds_from_genomic, fna_filepath, species, dir_out):
    global NUMBER_OF_NEED_TO_BE_WRITTEN
    global BROKEN_LIST
    try:
        ortho_data = pd.read_csv(orthodata_filepath, sep='\t', usecols=range(species + 1))  # +1 - 'Gene name' column
    except ValueError:
        ortho_data = pd.read_csv(orthodata_filepath, sep=',', usecols=range(species + 1))
    # for column_number, column_name in enumerate(ortho_data.columns):
    for column_name in ortho_data.columns:
        NUMBER_OF_NEED_TO_BE_WRITTEN = 0
        # species_numerating = str(column_number + 1)
        if column_name.isalnum():
            column_number = int(column_name)
            annotation_gbff_path = os.path.join(annotation_gbff, "{}.{}".format(column_name, 'gbff'))
            cds_from_genomic_path = os.path.join(cds_from_genomic, "{}.{}".format(column_name, 'fna'))
            logging.info('Working with species number {}'.format(column_name))
            """ define broken_list as a global var for re-extracting (post-processing) seqs from list ['item1', 'item2'..],
            e.g ['3557', '5781', '1503'] """
            # if BROKEN_LIST:
            #     broken_list_int_convert = list()
            #     for number in BROKEN_LIST:
            #         broken_list_int_convert.append(int(number) - 1)
            #         ortho_protein_ids = ortho_data.iloc[broken_list_int_convert, column_number].values
            # else:  # TODO: redo for gene_name instead of index
            ortho_protein_ids = ortho_data.loc[:, column_name].values
            gene_names = ortho_data.loc[:, 'Gene name'].values

            NUMBER_OF_NEED_TO_BE_WRITTEN = len(ortho_protein_ids)
            logging.info("NUMBER_OF_NEED_TO_BE_WRITTEN for species {} = {}".format(column_name,
                                                                                   NUMBER_OF_NEED_TO_BE_WRITTEN))
            get_and_write_nucleotide_seq(annotation_gbff_path, cds_from_genomic_path, ortho_protein_ids, gene_names,
                                         dir_out, column_name, fna_filepath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ortho', help='Path to the orthologs table file .tsv, see example in GScour/data,\'Gene '
                                        'name\' column is required', nargs='?', required=True)
    parser.add_argument('--gbff', help='Path to the folder with annotation .gbff files from '  # TODO: test necessity
                                       'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--cds', help='Path to the folder with refseq annotated files !.fna!'
                                      '_cds_from_genomic.fna or other source fasta format of'
                                      'the nucleotide sequences corresponding to all CDS'
                                      'features annotated on the assembly', nargs='?')
    # parser.add_argument('--csv', help='Path to the folder with annotation .csv files from '
    #                                  'www.ncbi.nlm.nih.gov/genome/', nargs='?')
    parser.add_argument('--genome', help='Path to the folder with .fna files', nargs='?')
    parser.add_argument('--species', help='Number of species', nargs='?', required=True)
    parser.add_argument('--group', help='Minimal size of species group', nargs='?', required=True)
    parser.add_argument('--out', help='Path to the folder for result write out, if it does not exist,'
                                      ' it will be created automatically', nargs='?', required=True)
    args = parser.parse_args()
    group = int(args.group)
    resulting_file = os.path.join(args.out, "get_ortho_nuc_result.xlsx")
    try:
        directory_out = args.out
        if not os.path.isdir(directory_out):
            os.makedirs(directory_out)
        main(args.ortho, args.gbff, args.cds, args.genome, int(args.species), directory_out)
        written_files_number = len(CORRECT_FILES)
        delta = NUMBER_OF_NEED_TO_BE_WRITTEN - written_files_number
        if delta == 0:
            logging.info("All files are written")

        for file, species_info in CORRECT_FILES.items():
            print("processed", file, species_info, group)
            if species_info[1] != group:
                if file not in BROKEN_SPECIES:
                    BROKEN_SPECIES.append(file)  # gene name
            else:
                "sort by folders-groups"
                CORRECT_RECORD.append(file)
                list_of_species = list(species_info[0].split(','))
                list_of_species.sort()
                group_folder_name = "".join(list_of_species)
                group_folder_path = os.path.join(directory_out, group_folder_name)
                if not os.path.isdir(group_folder_path):
                    os.makedirs(group_folder_path)
                os.replace(os.path.join(directory_out, file + ".fna"), os.path.join(group_folder_path, file + ".fna"))
        # resulting_file = os.path.join(args.out, "get_ortho_nuc_result.xlsx")
        writer = pd.ExcelWriter(resulting_file, engine='openpyxl')
        df_corr = pd.DataFrame({'Gene name': CORRECT_RECORD})
        df_corr.to_excel(writer, sheet_name='correct records', index=False)

        if BROKEN_SPECIES:
            broken_species_folder = os.path.join(args.out, "broken_species_files")
            os.makedirs(broken_species_folder)
            replace_broken_files_and_write_table(args.out, broken_species_folder, writer)
            residue = written_files_number - len(BROKEN_SPECIES)
            logging.info("NUMBER_OF_NEED_TO_BE_WRITTEN = {},  WRITTEN_FILES = {}, where {} in BROKEN_SPECIES list"
                         .format(NUMBER_OF_NEED_TO_BE_WRITTEN, written_files_number, len(BROKEN_SPECIES)))
            logging.info("removed broken species files into folder 'broken_species_files' in cwd,"
                         "please check out folder for .fna files number: {}".format(residue))
        writer.save()
    except BaseException as e:
        logging.exception("Unexpected error: {}".format(e))
    for k, v in ABSENT_IN_CDS.items():
        logging.warning("ABSENT_IN_CDS for {} species: {}, namely\n{}".format(k, len(ABSENT_IN_CDS[k]), v))
    logging.warning("BROKEN_SPECIES {} : {}".format(len(BROKEN_SPECIES), BROKEN_SPECIES))
    logging.warning("BROKEN_STOP_CODON {} : {}".format(len(BROKEN_STOP_CODON), BROKEN_STOP_CODON))
    logging.warning("BROKEN_START_CODON {} : {}".format(len(BROKEN_START_CODON), BROKEN_START_CODON))
    logging.warning("BROKEN_ACCORDANCE {} : {}".format(len(BROKEN_ACCORDANCE), BROKEN_ACCORDANCE))
    logging.warning("BROKEN_MULTIPLE_THREE {} : {}".format(len(BROKEN_MULTIPLE_THREE), BROKEN_MULTIPLE_THREE))
    logging.info("Correctly recorded files {}\n".format(len(CORRECT_RECORD)))
    logging.info("The work has been completed\nthe summary is in {}\nNucleotide sequences are in {}\n"
                 "".format(resulting_file, args.out))
    if BROKEN_SPECIES:
        logging.info("Broken nucleotide sequences are in {}".format(broken_species_folder))
