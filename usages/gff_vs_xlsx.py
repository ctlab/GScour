#!/usr/bin/env python
# -*- coding: utf-8 -*-
from BCBio import GFF
import pandas as pd

""" script to find intersection of some target genes (from .xlsx format) and genes from some species (annotation 
    .gff format)"""

in_file_gff = "/path_to/species.gff"
in_file_xlsx = "/path_to/genes_to_find.xlsx"
gff_genes = list()
in_handle = open(in_file_gff)
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        if feature.qualifiers.get('gene'):
            gene = feature.qualifiers.get('gene')[0]
            gff_genes.append(gene)
in_handle.close()
gff_genes_pd = set(gff_genes)
df = pd.io.excel.read_excel(in_file_xlsx, engine='openpyxl')
genes_immune = set(df['Gene symbol'])
intersection = pd.Series(list(gff_genes_pd.intersection(genes_immune)))

writer = pd.ExcelWriter("/home/alina_grf/progprojects/data/intersection.xlsx", engine='openpyxl')
intersection.to_excel(writer, sheet_name='intersection')
writer.save()
