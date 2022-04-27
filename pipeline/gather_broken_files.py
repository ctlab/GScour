import re
import pandas as pd

path_to_log = '/home/alina_grf/progprojects/GScour/pipeline/gblocks_alignment.log'
broken_file_path = '/home/alina_grf/progprojects/GScour/pipeline/broken_files_gblocks.xlsx'
pattern = re.compile('ERROR : Infile (\d+)')
broken_files = list()

with open(path_to_log, 'r') as f:
    for line in f:
        if re.search(pattern, line):
            broken_file_number = re.search(pattern, line).group(1)
            broken_files.append(broken_file_number)


writer = pd.ExcelWriter(broken_file_path, engine='xlsxwriter')
df = pd.DataFrame({'broken files': broken_files})
df.to_excel(writer)
writer.save()



