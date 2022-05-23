import argparse
import pandas as pd

"""
script co compare data from previous research with current data, all tables must have column 'Gene name'
"""
research_sheet_names = ['Datasheet S6', 'Datasheet S5a']


def main(data_sheet_path, adjust_sheet_path, out_sheet_path):
    adjust_sheet = pd.io.excel.read_excel(adjust_sheet_path, sheet_name=None, engine='openpyxl')
    data_sheet = pd.io.excel.read_excel(data_sheet_path, sheet_name=None, engine='openpyxl')
    writer = pd.ExcelWriter(out_sheet_path, engine='openpyxl')
    for name_data, sheet_data in data_sheet.items():
        if name_data in research_sheet_names:
            full_out_table = pd.DataFrame()
            for name_adjust, sheet_adjust in adjust_sheet.items():
                merged = pd.merge(sheet_data, sheet_adjust, how='inner', on=['Gene name'],
                                  suffixes=("_article", "_new"), indicator=True)
                d = {
                    "left_only": "Only_article", "right_only": "Only_new",
                    "both": "{}".format(name_adjust)
                    }
                merged['_merge'] = merged['_merge'].map(d)
                full_out_table = full_out_table.append(merged)
            if name_data == 'Datasheet S6':
                full_out_table.to_excel(writer, sheet_name='Datasheet S6', index=False)
            elif name_data == 'Datasheet S5a':
                full_out_table.to_excel(writer, sheet_name='Datasheet S5a', index=False)
    writer.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', help='The path to the data .xlsx to be compared with', nargs='?', required=True)
    parser.add_argument('--i', help='Path to the .xlsx file with sheets with column \'Gene name\'',
                        nargs='?', required=True)
    parser.add_argument('--out', help='Path to the output .xlsx file with intersection of data', nargs='?',
                        required=True)
    args = parser.parse_args()
    try:
        main(args.data, args.adjust, args.out)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
    print("The work has been completed")
