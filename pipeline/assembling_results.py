#!/usr/bin/env python
import argparse
import logging
import os
import re

import pandas as pd

""" gathering results from common_sheets, remove duplications"""

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def parse_dir(in_directory):
    for f in os.scandir(in_directory):
        if f.name.split('.')[-1] == 'xlsx':
            yield f.name


def main(input_folder):
    summary_sheets = []
    for table in parse_dir(input_folder):
        logger.info("Work with table {}".format(table))
        summary_sheet = pd.io.excel.read_excel(os.path.join(input_folder, table), sheet_name='summary', index_col=0,
                                               engine='openpyxl')
        summary_sheets.append(summary_sheet)
    common_summary = pd.concat(summary_sheets)
    common_summary.drop_duplicates(subset=['Gene name', 'p-value', 'NCBI protein_id', 'Species group'], inplace=True,
                                   ignore_index=True)

    for gene_name in common_summary['Gene name']:
        try:
            pure_gene_name = re.search(r'([a-zA-Z0-9]+)_', gene_name).group(1)
            common_summary['Gene name'].replace(to_replace=gene_name, value=pure_gene_name, inplace=True)
            print(common_summary['Gene name'])
        except:
            print('error with ', gene_name)
            continue
    """
    there are variants:
    common_summary.drop_duplicates(subset=['Gene name', 'Species group'], inplace=True,
                                   ignore_index=True)
    common_summary.drop_duplicates(subset=['Gene name', 'NCBI protein_id', 'Species group'], inplace=True,
                                   ignore_index=True)
    """
    common_summary_aggregation = common_summary.groupby(['Gene name']).agg({'p-value': list, 'NCBI protein_id': list,
                                                                           'Species group': list}).reset_index()

    common_summary_min = common_summary.groupby(['Gene name', 'NCBI protein_id', 'Species group'])['p-value'].min()
    # common_summary['Species group'] = common_summary['Species group'].astype(str)
    # common_summary = common_summary.groupby('Gene name').agg({'Species group': ', '.join}).reset_index()
    writer = pd.ExcelWriter(os.path.join(input_folder, '{}.xlsx'.format("assembled_results")), engine='openpyxl')
    common_summary_aggregation.to_excel(writer, sheet_name='from all groups')
    common_summary_min.to_excel(writer, sheet_name='chosen min p-val')

    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', help='Path to the folder with final common sheets with sheet \'summary\'',
                        nargs='?', required=True)
    args = parser.parse_args()
    logger.info("Passed args: input directory with data {}".format(args.data))
    try:
        main(args.data)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
    logger.info("The work has been completed")