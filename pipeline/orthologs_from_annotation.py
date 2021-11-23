#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from BCBio import GFF
import pandas as pd

""" script to find orthologs among annotation .gff files for some set of species """

annotation_path_folder = '/path/to/gff_files'
result_xlsx_file_path = "/path/to/orthologs.xlsx"
gff_genes = dict()
species_names = list()
species_record = dict()


def parse_folder(dir_path):
    for f in os.scandir(dir_path):
        if f.name.split('.')[-1] == 'gff':
            yield f, f.name.split('.')[0]


def get_genes_from_gff():
    for gff_file, species_name in parse_folder(annotation_path_folder):
        species_record[species_name] = list()
        species_names.append(species_name)
        in_handle = open(gff_file)
        gff_genes[species_name] = set()
        for rec in GFF.parse(in_handle):
            rec_id = rec[0]
            print("rec_id=".format(rec_id))
            for feature in rec.features:
                if feature.qualifiers.get('gene'):
                    gene = feature.qualifiers.get('gene')[0]
                    gff_genes[species_name].add(gene)
                    species_record[species_name].append((gene, rec_id))
        print("for {}: {} genes".format(species_name, len(gff_genes[species_name])))
        in_handle.close()


def pairwise_comparison():
    for i in range(1, len(species_names)):
        print("Error prevention - pairwise comparison: intersection of {} and {} = {}".format(species_names[0],
                                                                                              species_names[i],
                                                                                              len(set.intersection(
                                                                                                  gff_genes[
                                                                                                      species_names[0]],
                                                                                                  gff_genes[
                                                                                                      species_names[
                                                                                                          i]]))))


def find_orthologs():
    sets = [v for v in gff_genes.values()]
    orthologs_intersection = set.intersection(*sets)
    return orthologs_intersection


def write_result(result, result_file_path):
    result = pd.Series(list(result))
    writer = pd.ExcelWriter(result_file_path, engine='openpyxl')
    result.to_excel(writer, sheet_name='orthologs')
    writer.save()


if __name__ == '__main__':
    get_genes_from_gff()
    pairwise_comparison()
    orthologs = find_orthologs()
    write_result(orthologs, result_xlsx_file_path)
    print("done")
