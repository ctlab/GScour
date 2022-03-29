#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from BCBio import GFF
import pandas as pd
from Bio import SeqIO

""" script to find orthologs among annotation .gff or .gbff files for some set of species """

annotation_path_folder = '/path/to/gff_files'
result_xlsx_file_path = "/path/to/orthologs.xlsx"
genes = dict()
species_names = list()


def parse_folder(dir_path, annotation_type):
    for f in os.scandir(dir_path):
        if annotation_type == 'gff':
            if f.name.split('.')[-1] == 'gff':
                yield f, f.name.split('.')[0]
        elif annotation_type == 'gbff':
            if f.name.split('.')[-1] == 'gbff':
                yield f, f.name.split('.')[0]


def get_genes_from_gff(annotation_file_path, species_name):
    in_handle = open(annotation_file_path)
    for rec in GFF.parse(in_handle):
        rec_id = rec[0]
        print("rec_id=".format(rec_id))
        for feature in rec.features:
            if feature.qualifiers.get('gene'):
                gene = feature.qualifiers.get('gene')[0]
                genes[species_name].add(gene)
    print("for {}: {} genes".format(species_name, len(genes[species_name])))
    in_handle.close()


def get_genes_proteins_from_gbff(annotation_file_path, species_name):
    gene_protein_len = dict()
    for record in SeqIO.parse(annotation_file_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers.get("protein_id"):
                    protein_id = feature.qualifiers.get("protein_id")[0]  # str
                    protein_translation = feature.qualifiers.get("translation")[0]  # str
                    protein_translation_length = len(protein_translation)  # int
                    gene = feature.qualifiers.get('gene')[0]
                    if not gene_protein_len.get(gene):
                        gene_protein_len[gene] = {"protein_id": protein_id, "length": protein_translation_length}
                    else:
                        if gene_protein_len[gene]["length"] < protein_translation_length:
                            print("the current protein length is less than new, replace: old {}".format(
                                gene_protein_len[gene]))
                            gene_protein_len[gene] = {"protein_id": protein_id, "length": protein_translation_length}
                            print("for new {}".format(gene_protein_len[gene]))
    for k, v in gene_protein_len.items():
        print("k, v = {}, {}".format(k, v))
        genes[species_name].add(k)
    print("for species number {}: {} genes".format(species_name, len(genes[species_name])))
    return gene_protein_len


def get_genes_from_annotation(annotation_type):
    gene_protein_for_species = dict()
    for annotation_file_path, species_name in parse_folder(annotation_path_folder, annotation_type):
        species_names.append(species_name)
        genes[species_name] = set()
        if annotation_type == 'gff':
            get_genes_from_gff(annotation_file_path, species_name)
        if annotation_type == 'gbff':
            gene_protein_len = get_genes_proteins_from_gbff(annotation_file_path, species_name)
            gene_protein_for_species[species_name] = gene_protein_len
    return gene_protein_for_species


def pairwise_comparison():
    for i in range(1, len(species_names)):
        print("Smoke test, comparison: intersection of {} and {} = {}".
              format(species_names[0], species_names[i], len(set.intersection(genes[species_names[0]],
                                                                              genes[species_names[i]]))))


def extract_orthologs():
    gene_vals = [v for v in genes.values()]
    orthologs_intersection = set.intersection(*gene_vals)
    print("orthologs intersection of length {}".format(len(orthologs_intersection)))
    return orthologs_intersection


def write_result(result, result_file_path, gene_protein_dict):
    if gene_protein_dict:
        gene_names = dict()
        proteins = dict()
        species_names = list()
        for species, vals in gene_protein_dict.items():
            species_names.append(species)
            gene_names[species] = list()
            proteins[species] = list()
            for gene in result:
                if vals.get(gene):
                    gene_names[species].append(gene)
                    proteins[species].append(vals[gene]['protein_id'])

        # header = pd.MultiIndex.from_product([species_names,
        #                                      ['gene', 'protein_id']])
        df1 = pd.DataFrame(gene_names)
        df2 = pd.DataFrame(proteins)
        dfs = [df1, df2]
        res = pd.concat(dfs, axis=1)
        result = pd.DataFrame(res)
    else:
        result = pd.Series(list(result))
    writer = pd.ExcelWriter(result_file_path, engine='openpyxl')
    result.to_excel(writer, sheet_name='orthologs', index=False)
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', help='Format of annotation files gbff or gff', choices=['gbff', 'gff'],
                        required=True, nargs='?')
    args = parser.parse_args()
    gene_protein_for_species = get_genes_from_annotation(args.f)
    pairwise_comparison()
    orthologs = extract_orthologs()
    write_result(orthologs, result_xlsx_file_path, gene_protein_for_species)
    print("done")
