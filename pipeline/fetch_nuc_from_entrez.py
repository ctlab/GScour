import os
from Bio import Entrez
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

out_fasta_path = os.path.join(os.getcwd(), 'fetched.fasta')
out_fasta_log = os.path.join(os.getcwd(), 'fetched.log')
print("out_file_path", out_fasta_path)
Entrez.email = "A.N.Other@example.com"
"""
script to fetch nucleotide sequences from ncbi nucleotide db by NCBI Reference Sequence IDs,
for example:
"""

ids_to_fetch = ['NM_001300876.1', 'XM_004752454.2', 'XM_007082235.2', 'XM_015063672.2', 'XM_029955397.1',
                'XM_039245653.1', 'XM_004405986.1', 'XM_032347068.1', 'XM_032854848.1', 'XM_026500060.1',
                'XM_019455066.1', 'XM_023257226.1', 'XM_025926639.1', 'XM_040491460.1', 'XM_030323191.1',
                'XM_022494882.1', 'XM_032420589.1', 'NM_001290067.1', 'XM_036110399.1', 'XM_021678605.1',
                'XM_035002168.1', 'XM_027594348.2', 'XM_028110796.1', 'XM_025882721.1', 'NM_001394523.1',
                'XM_008684278.2', 'XM_026001743.1', 'XM_041754114.1', 'XM_025476634.2']
gene_name = 'LYZ'
print("ids_to_fetch:", len(ids_to_fetch))


def rename_seq(name):
    if re.search(r'\s([a-zA-Z, ]+)\s\([A-Za-z0-9]+\)', name):
        found = re.search(r'\s([a-zA-Z, ]+)\s\([A-Za-z0-9]+\)', name).group(1)
        res = found.replace(' ', '_')
        return res
    else:
        print('not matched name', name)


def fetch_fasta():
    for nuc in ids_to_fetch:
        handle = Entrez.efetch(db="nucleotide", id=nuc, rettype="fasta", retmode="text")
        result = handle.read()
        handle.close()
        yield result, nuc


def main():
    species_counter = 0
    species_dict = dict()
    with open(out_fasta_log, 'a') as log_file:
        log_header = 'gene_name id seq_length species_numerating\n'
        log_file.write(log_header)
        for string_seq, seq_id in fetch_fasta():
            print("string_seq\n", string_seq)
            species_counter += 1
            splitted = string_seq.split('\n')
            species_name = rename_seq(splitted[0])
            species_dict[species_name] = str(species_counter)
            nuc_seq = ''.join(splitted[1:])
            seq_length = len(nuc_seq)
            seq_record = SeqRecord(Seq(nuc_seq), id=str(species_counter), description='')
            print(seq_record)
            with open(out_fasta_path, "a") as f:
                SeqIO.write(seq_record, f, "fasta")
                log_file.write('{} {} {} {}\n'.format(gene_name, seq_id, seq_length, species_counter))
        print("species counter=", species_counter, "assigned number to species: ", repr(species_dict))


main()
