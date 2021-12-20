import gseapy as gp
import pandas as pd

in_table_path = ''
out_table_path = ''
summary_sheet = pd.io.excel.read_excel(in_table_path, sheet_name='summary', index_col=0,
                                               engine='openpyxl')
gene_names = summary_sheet['Gene name']
gene_list = gene_names.squeeze().str.strip().tolist()
print(gene_list[:10])
libraries = gp.get_library_name()  # default: Human
enr = gp.enrichr(gene_list=gene_list, description='pathway', gene_sets=libraries, outdir='test')
#
print(type(enr.results.head(5)), enr.results.head(5))

writer = pd.ExcelWriter(out_table_path, engine='xlsxwriter')
enr.results.to_excel(writer, sheet_name='GSEA')
writer.save()