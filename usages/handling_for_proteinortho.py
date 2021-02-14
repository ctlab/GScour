import re

"""
small script for corrected raw genome data (.gff) for processing Proteinortho - 
added appropriate Name and next writing

from specs for POFF Proteinortho:
As Proteinortho is primarily made for proteins, it will only accept GFF
entries of type CDS (column #3 in the GFF-file). The attributes column
(#9) must contain Name=GENE IDENTIFIER where GENE IDENTIFIER corresponds
to the respective identifier in the FASTA format.

+ .gff file (#9) must contain Name=GENE IDENTIFIER then a semicolon (;) and next writing,
for example ";Note=ExampleNote" to be fed by Proteinortho 
"""

in_file = "/example_home/example_in.gff"
out_file = "/example_home/example_out.gff"

stencil_name = r'([a-z, _]+-[A-Z, \d, \.]+-[a-z]+-[a-z]+-\d+\.\d+-[a-z, A-Z]+-\d+)'  # names for raw lion genome

with open(in_file, 'r') as in_f:
    with open(out_file, 'w+', encoding='utf-8') as o_f:
        for line in in_f:
            if re.search('CDS', line) and re.search(stencil_name, line):
                lines = line.split('\t')
                name = re.search(stencil_name, line).group(0)
                new_line = lines[8][:-1]
                new_line += ';Name={};Note=ExampleNote\n'.format(name)
                lines[8] = new_line
                for li in lines[:-1]:
                    o_f.write(li + '\t')
                o_f.write(lines[-1])
            else:
                o_f.write(line)