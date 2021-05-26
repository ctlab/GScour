import argparse
import pandas as pd


def get_bonferroni_significance(row):
    if row['P-value'] < row['Bonferroni']:
        return '+'
    else:
        return '-'


def get_FDR_significance(row):
    if row['P-value'] < row['FDR']:
        return '+'
    else:
        return '-'


def main(sheet_path, out_path):
    common_sheet = pd.io.excel.read_excel(sheet_path, sheet_name=None, index_col=0, engine='openpyxl')
    writer = pd.ExcelWriter(out_path, engine='openpyxl')
    for name, sheet in common_sheet.items():
        sheet.sort_values(by="P-value", ignore_index=True, inplace=True)
        sheet['Bonferroni'] = 0.05 / sheet["P-value"].shape[0]
        sheet['Significance_Bonf'] = sheet.apply(lambda row: get_bonferroni_significance(row), axis=1)
        sheet['FDR'] = (sheet["P-value"].index + 1) * 0.05 / sheet["P-value"].shape[0]
        sheet['Significance_FDR'] = sheet.apply(lambda row: get_FDR_significance(row), axis=1)
        sheet.to_excel(writer, sheet_name=name)
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sheet', help='Path to the common_sheet.xlsx (with removed last sheet \'summary\')'
                                        ' - result of \'paml_out_analysis.py\'', nargs='?', default="codeml")
    parser.add_argument('--out', help='Path to the output.xlsx file', nargs='?', default="codeml")
    args = parser.parse_args()
    try:
        main(args.sheet, args.out)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
    print("The work has been completed")