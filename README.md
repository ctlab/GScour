# Genome-scale detection of positive selection
## Requirements
### Perl
Perl 5: https://www.perl.org/
### Aligners and alignment analysis tools
PRANK multiple sequence aligner (http://wasabiapp.org/software/prank/) - v.140603
GUIDANCE (http://guidance.tau.ac.il/) - v1.5
Note that a bug fix is required for GUIDANCE (version 1.5 - 2014, August 7) to work with PRANK, see GUIDANCE_source_code_fix_for_running_PRANK
PAML software package, which includes codeml (http://abacus.gene.ucl.ac.uk/software/paml.html) - v4.8a

## Steps
### 1. One-to-one orthologs
Obtain one-to-one ortholog clusters for whole-genome sequences.

### 2. Sequences
#### 2a. Get sequences

### 3. Alignments
Produce codon-based nucleotide sequence alignments for all the one-to-one ortholog clusters. Then assess the confidence in the alignments using two independent approaches.

#### 3a. PRANK codon-based multiple alignment
#### 3b. GUIDANCE assessment and masking
NOTE: this step takes a lot of computation time.

### 4. Evolutionary analyses
Perform maximum likelihood (ML) dN/dS analysis to infer positive selection of genes and codons, using codeml from the PAML software package.
