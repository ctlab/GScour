# Genome-scale detection of positive selection
## Requirements
### Perl
Perl 5: https://www.perl.org/
### Aligners and alignment analysis tools
##### PRANK
multiple sequence aligner (http://wasabiapp.org/software/prank/) - v.170427
For launch scripts without significant changes add path to prank location in your environment variable PATH like
export PATH="$PATH:/path/to/dir/prank/bin/"

GUIDANCE (http://guidance.tau.ac.il/) - v2.02
GUIDANCE allows using alignment (MAFFT, PRANK, CLUSTALW) as a subprocess.

## Steps
### 1. One-to-one orthologs
Obtain one-to-one ortholog clusters for whole-genome sequences.

### 2. Sequences
#### 2a. Get sequences

### 3. Alignments
Produce codon-based nucleotide sequence alignments for all the one-to-one ortholog clusters. Then assess the confidence in the alignments using two independent approaches.

#### 3a. PRANK codon-based multiple alignment
PRANK may be used as a subprocess of GUIDANCE, therefore skip this step.
#### 3b. GUIDANCE assessment and masking
NOTE: this step takes a lot of computation time.
If you use alignment (MAFFT, PRANK, CLUSTALW) as a subproces of GUIDANCE you should comment out the lines 
launch = 'perl guidance.pl --seqFile {0} --msaProgram PRANK --seqType nuc --msaFile {0} --proc_num {1} --outDir {2}'.format(infile, threads, outdir)
in run_guidance_alignment.py and set the desired value instead of PRANK.


### 4. Evolutionary analyses
Perform maximum likelihood (ML) dN/dS analysis to infer positive selection of genes and codons, using codeml from the PAML software package.
