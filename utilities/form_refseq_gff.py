#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import re
import os

"""
extract from gff file
protein-coding genes and CDS types
obtained by methods for RefSeq  
"""

refseq_gene = [r'\tRefSeq\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding',
               r'\tBestRefSeq\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding',
               r'\tCurated Genomic\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding',
               r'\tGnomon\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding',
               r'\tBestRefSeq%2CGnomon\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding',
               r'\tCurated Genomic%2CGnomon\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding',
               r'\ttRNAscan-SE\tgene\t\d+\t\d+\t(.*?)\t(.*?)\t(.*?)\t(.*?)gene_biotype=protein_coding']
refseq_cds = [r'\tRefSeq\tCDS\t',
              r'\tBestRefSeq\tCDS\t',
              r'\tCurated Genomic\tCDS\t',
              r'\tGnomon\tCDS\t',
              r'\tBestRefSeq%2CGnomon\tCDS\t',
              r'\tCurated Genomic%2CGnomon\tCDS\t',
              r'\ttRNAscan-SE\tCDS\t']


def main(initial_dir, corrected_dir):
    for gff_file in os.listdir(initial_dir):
        if gff_file.endswith(".gff"):
            counter_gene = 0
            counter_cds = 0
            with open(os.path.join(initial_dir, gff_file), 'r') as in_f:
                with open(os.path.join(corrected_dir, gff_file), 'w+') as out_f:
                    for line in in_f:
                        if re.search('^#', line):
                            out_f.write(line)
                            continue
                        for refseq_method, cds in zip(refseq_gene, refseq_cds):
                            if re.search(refseq_method, line):
                                counter_gene += 1
                                out_f.write(line)
                                break
                            if re.search(cds, line):
                                counter_cds += 1
                                out_f.write(line)
                                break

            print("{} has been wrote in corrected dir, counter_gene = {}, counter_cds = {}".
                  format(gff_file, counter_gene, counter_cds))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='The full path to the folder with .gff files', nargs='?', required=True)
    args = parser.parse_args()
    gff_initial_dir = args.i
    gff_corrected_dir = os.path.join(gff_initial_dir, "gff_corr")
    if not os.path.isdir(gff_corrected_dir):
        os.makedirs(gff_corrected_dir)
        print("{} directory has been created for corrected .gff files".format(gff_corrected_dir))
    try:
        main(gff_initial_dir, gff_corrected_dir)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
