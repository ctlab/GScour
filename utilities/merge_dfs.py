import pandas as pd

"""small script ro merge or concat dataframes, see pandas docs for details"""

path1 = '/path/formed_acfhp.poff.tsv'
path2 = '/path/orthologs_from_annotation.xlsx'
result_file_path = '/abspath/final_table_concat_poff_and_annotation.tsv'

df1 = pd.read_csv(path1, sep='\t')
df2 = pd.read_excel(path2, engine='openpyxl')

# res = df1.merge(df2, on='Gene name', how='outer')
res = pd.concat([df1, df2])
print("number of samples before removing duplications={}".format(res.shape[0]))
res.drop_duplicates(subset=['1', '2', '3', '4', '5'], inplace=True)
print("number of samples after removing duplications={}".format(res.shape[0]))
res.to_csv(result_file_path, sep="\t")
