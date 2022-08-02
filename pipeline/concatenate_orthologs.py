import pandas as pd

"""small script ro merge or concat dataframes, see pandas docs for details"""

path1 = 'formed_acfhp.poff.tsv'
path2 = 'orthologs_from_annotation.xlsx'
result_file_path = 'final_orthologs_table_concatenation.tsv'

df1 = pd.read_csv(path1, sep='\t')
df2 = pd.read_excel(path2, engine='openpyxl')


def rename_duplications(df, suffix='-duplicate-'):
    df.loc[df.duplicated('Gene name', keep=False), 'Gene name'] = df.loc[df.duplicated('Gene name',
                                                                                       keep=False), 'Gene name'] + \
                                                                  suffix + (df[df.duplicated('Gene name',
                                                                                             keep=False)].groupby(
                                                                            "Gene name").cumcount() + 1).astype(str)
    return df


# res = df1.merge(df2, on='Gene name', how='outer')
res = pd.concat([df1, df2])
print("number of samples before removing duplications={}".format(res.shape[0]))
res.drop_duplicates(subset=['1', '2', '3', '4', '5'], inplace=True)  # remove duplications by protein samples
print("number of samples after removing duplications of protein samples ={}".format(res.shape[0]))
# res.set_index('Gene name', inplace=True)
# print("idx", res.index)
res = rename_duplications(res)  # rename duplications by 'Gene name'
res.to_csv(result_file_path, sep="\t")
