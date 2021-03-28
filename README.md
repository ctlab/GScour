# Genome-scale detection of positive selection
## Requirements
- Biopython
- PAML (Phylogenetic Analysis by Maximum Likelihood)
- Numpy
- Pandas
- Proteinortho (https://gitlab.com/paulklemm_PHD/proteinortho)
### Aligners and alignment analysis tools
- PRANK
Multiple sequence aligner (http://wasabiapp.org/software/prank/) - v.170427
For launch scripts without significant changes add path to prank location in your environment variable PATH like
export PATH="$PATH:/path/to/dir/prank/bin/"
- GUIDANCE (http://guidance.tau.ac.il/) - v2.02
GUIDANCE allows using alignment (MAFFT, PRANK, CLUSTALW) as a subprocess.
- Gblocks (http://gensoft.pasteur.fr/docs/gblocks/0.91b/)
- SWAMP (https://github.com/peterwharrison/SWAMP)
## Steps
### 1. One-to-one orthologs
Obtain one-to-one ortholog clusters for whole-genome sequences.
#### 1.1 Finding orthologous proteins
- It is necessary for futher analysis to name the species in numbers (1,2,3,...);
- Launch proteinortho. For example, launch_proteinortho.sh or launch_proteinortho_synteny.sh.
Result file will be needed: project_name.proteinortho.tsv or project_name.poff.tsv if synteny;
#### 1.2 Form ortologs table
- Form a ortologs table with fulfillment of requirements: single-copy orthologs only, group (minimal
number of species in group), required species (one target species in relation to which the analysis is made).
Use --help for help with arguments. Result will be recorded to project_name_formed_orthologs_table.tsv;  
`python get_orthologs_table.py --project project_name --species 8 --group 6 --required 6`

### 2. Sequences
#### 2.1. Get nucleotide sequences
Annotation .gbff, genomes .fna and cds_from_genomic.fna are needed. Use --help for help with arguments.
Firstly, try to extract sequences from cds_from_genomic.fna, else from .gbff. Result: .fna file with sequences and
.log file with summary for every .fna. See log in "get_ortho_nuc_seqs.log".  
`python get_ortho_nucleotides.py --ortho /abspath/to/thetable/project_name_formed_orthologs_table.tsv --gbff /abspath/tothe/gbff_folder/gbff/ --cds /abspath/tothe/cds/cds_refseq/ --genome /abspath/tothe/fnafiles/genomes/ --species 8 --group 6 --out /abspath/tothe/nuc_out_folder/`
#### 2.2 Check duplicates
Perform additional check to exclude duplicates  
`python check_duplicates.py --i /abspath/tothe/nuc_out_folder(from_step2.1)/`

### 3. Alignments
Produce codon-based nucleotide sequence alignments for all the one-to-one ortholog clusters.
#### 3.1 PRANK codon-based multiple alignment
PRANK may be used separetly (recommended) or as a subprocess of GUIDANCE (therefore skip this step).
See --help for arguments. Option --tree isn't adapted to work with groups, therefore it should be used if there is only
one group (args in get_orthologs_table.py group==species).  
`python prank_alignment.py --i /abspath/tothe/nuc_out_folder --o /abspath/tothe/nuc_out_prank/ --threads 32`
#### 3.2. GUIDANCE assessment and masking
NOTE: this step takes a lot of computation time.
If you use MAFFT or CLUSTAL (not PANK) as a subproces of GUIDANCE you should correct --msaProgram option
in the line https://github.com/ctlab/search_for_positive_selection/blob/4748195803e635284e77007375e2b699db922cbb/pipeline/guidance_alignment.py#L47  
See --help for help with arguments.  
`python guidance_alignment.py --i /abspath/tothefolder/with_nucseqs/ --o /abspath/tothe/guidance_out/ --exec /abdpath/guidance.v2.02/www/Guidance/guidance.pl --threads 22`


### 4. Evolutionary analyses
Perform maximum likelihood (ML) dN/dS analysis to infer positive selection of genes and codons, using codeml from the PAML software package.
