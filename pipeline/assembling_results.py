#!/usr/bin/env python
import argparse
import logging
import os
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
    common_summary.drop_duplicates(subset='Gene name', inplace=True, ignore_index=True)
    writer = pd.ExcelWriter(os.path.join(input_folder, '{}.xlsx'.format("assembled_results")), engine='openpyxl')
    common_summary.to_excel(writer, sheet_name='assembled')
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', help='Path to the folder with final common sheets with sheet \'summary\'',
                        nargs='?')
    args = parser.parse_args()
    logger.info("Passed args: input directory with data {}".format(args.data))
    try:
        main(args.data)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
    logger.info("The work has been completed")