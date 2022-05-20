#!/usr/bin/env python
import re
import pandas as pd

path_to_log = 'path'
broken_file_path = 'path'
pattern = re.compile('ERROR : Infile (\d+)')
broken_files = list()


with open(path_to_log, 'r') as f:
    for line in f:
        if re.search(pattern, line):
            broken_file_number = re.search(pattern, line).group(1)
            broken_files.append(broken_file_number)

if broken_files:
    writer = pd.ExcelWriter(broken_file_path, engine='xlsxwriter')
    df = pd.DataFrame({'broken files': broken_files})
    df.to_excel(writer)
    writer.save()
else:
    print("Luck. No broken files")



