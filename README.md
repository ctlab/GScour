# Genome-wide detection of positive selection
## Requirements
- Biopython
- Numpy
- Pandas
### Find orthologs
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
### Phylogenetic Analysis
- PAML (Phylogenetic Analysis by Maximum Likelihood)
## GScour flow
![GitHub Logo](https://github.com/ctlab/GScour/blob/master/GScour%20flow.jpg)
## Steps
See --help for help with arguments for any python script, logs are named as scripts with .log extension or just standard output. 
Output directory (with --o or --out option) will be created automatically and hold all output files [please provid full (and not relative) path]<br />
Errors can be found by keyword "WARNING" in logs.<br />
If there are any questions, errors, suggestions feel free to contact me via email
lnfyodorova@gmail.com.
### 1. One-to-one orthologs
Obtain one-to-one ortholog clusters for whole-genome sequences.
#### 1.1 Find orthologs proteins
- It is necessary for futher analysis to name the species in numbers (1,2,3,...) and to name all associated files in numbers  
(1.faa, 2.faa..., 1.gbff, 2.gbff...)
##### 1.1.1 Proteinortho
 - launch_proteinortho
For example, launch_proteinortho.sh or launch_proteinortho_synteny.sh.
Result file will be needed for futher analysis: project_name.proteinortho.tsv or project_name.poff.tsv if use synteny;
 - form ortologs table
Form a ortologs table with fulfillment of requirements: single-copy orthologs only, group (minimal
number of species in group), required species (one target species in relation to which the analysis is made). Result will be recorded to the file with prefix 'formed', for example, formed_project_name.poff.tsv;  
`python pipeline/get_orthologs_table.py --ortho /path/tothe/project_name.poff.tsv --species 8 --group 6`<br />` --required 6`
##### 1.1.2 Extract from ncbi
If all of needed species have ncbi annotations, we can extract orthologs from it.<br />
This requires defining variables in `pipeline/orthologs_from_annotation.py` script: <br />
- annotation_path_folder (folder with .gff files for every species)<br />
- result_file_path (path to the result .xlsx file that will be written)<br />
##### 1.1.3 Combine the previous ways
Merge or concatenate results, write paths and options into script and run
`python utilities/merge_dfs.py'
### 2. Sequences
#### 2.1. Get nucleotide sequences
##### 2.1.1 Extract sequences in accordance with orthologs result file
See example of file data/orthologs_table.tsv. 
Annotation .gbff or cds_from_genomic.fna are needed, genomes .fna are used for additional check. If you do not have one of these options, give path to empty folder. Result: .fna files with sequences and .log files with summary for every set of orthologs. See log in "get_ortho_nuc_seqs.log".  
`python pipeline/get_ortho_nucleotides.py --ortho /abspath/to/thetable/project_name_formed_orthologs_table.tsv`<br />`--gbff /abspath/tothe/gbff_folder/gbff/ --cds /abspath/tothe/cds/cds_refseq/`<br />`--genome /abspath/tothe/fnafiles/genomes/ --species 8 --group 6 --out /abspath/tothe/nuc_out_folder/`<br />
Sequences will be sorted to subfolders with unique names (group names)
corresponding to set of species in fasta file (sorted in increasing order). For example:  
```bash
$ cd /abspath/tothe/nuc_out_prank/  
$ ls  
1.fasta 2.fasta 3.fasta  
$ less 1.fasta    $ less 25.fasta    $ less 34.fasta  
>1                >2                >3  
ATG...            ATG...            ATG...  
>2                >3                >2      
ATG...            ATG...            ATG...
```

Result:  
```bash
$ cd /abspath/tothe/nuc_out_prank/  
$ ls */  
12/:            23/:  
1.fasta         25.fasta    
                34.fasta  
 ```                    
 ##### 2.1.2 Extract sequences in accordance with some target gene names list
`python pipeline/get_nucleotides_target_genes.py --t /path/to/orthologs.xlsx --gbff /path/to/folder_gbff_annotations` <br />
`--o /path/to/output_folder`
#### 2.2 Check duplicates
Perform additional check to exclude duplicates within one sample<br />
`python utilities/check_duplicates.py --i /abspath/tothe/nuc_out_folder(from_step2.1)/`

### 3. Alignments
Produce codon-based nucleotide sequence alignments for all the one-to-one ortholog clusters.
#### 3.1 PRANK multiple alignment
PRANK may be used separetly or as a subprocess of GUIDANCE.
Option --tree isn't adapted to work with groups, therefore it should be used if there is only one group (i.e args in get_orthologs_table.py group==species). By default we use option -translate for prank, but you can change it for -codon: codon alignment produces more accurate alignments than alignment of translated protein sequences. Whether you use tree or set output format this change can be made in 
https://github.com/ctlab/GScour/blob/5f1a4463ee29b0f4ef9f80cefc8d74c73e324868/pipeline/prank_alignment.py#L45
or lines 48, 51, 56.<br />
`python pipeline/prank_alignment.py --i /abspath/tothe/nuc_out_folder --o /abspath/tothe/nuc_out_prank/ --threads 32`
#### 3.2. GUIDANCE, masking of inconsistent residues
NOTE: this step can take a lot of computation time.
If you use MAFFT or CLUSTAL (not PANK) as a subproces of GUIDANCE you should correct --msaProgram option
in the line 
https://github.com/ctlab/GScour/blob/4748195803e635284e77007375e2b699db922cbb/pipeline/guidance_alignment.py#L47<br />
`python pipeline/guidance_alignment.py --i /abspath/tothefolder/with_nucseqs/ --o /abspath/tothe/guidance_out/` <br />
`--exec /abdpath/guidance.v2.02/www/Guidance/guidance.pl --threads 22`<br />
The resulting files stored in cleansed folder `/abspath/tothe/guidance_out/cleansed/`.
#### 3.3 Sort by groups
For example:<br />
`python utilities/sort_by_groups.py --i /abspath/tothe/nuc_out_prank/`<br />
Format of input file can be adjusted here <br />
https://github.com/ctlab/GScour/blob/a20f24a45a2a6163cbbb4834c726395d59438933/utilities/sort_by_groups.py#L48
#### 3.4 Gblocks, select conserved blocks of sequence
Use parameter --auto for automatic selection of gblocks parameters based on number of sequences for each group or adjust parameters to your needs in the params_string:  
https://github.com/ctlab/GScour/blob/7bd285734a26c521a844d08b8e4adcfa22804744/pipeline/gblocks_alignment.py#L35 <br />
`python pipeline/gblocks_alignment.py --i /abspath/tothe/nuc_out_prank/ --auto y --exec /abspath/Gblocks_0.91b/Gblocks`<br />`--threads 2`

### 4. Evolutionary analysis
#### 4.1 Provide trees
Phylogenetic trees must be provided for every species group in format as required for PAML.  
Name of tree should be the same as name of species folder ('12' -> '12.tree'). Put trees to separate folder. 
#### 4.2 Preprocessing, set right order for paml
Step can be skipped to 4.3.1 if the order is known.
##### 4.2.1 Replace files for test order
`python utilities/replace_for_test_order.py --i /abspath/tothe/nuc_out_prank/ --o /abspath/tothe/test_order/` 
##### 4.2.2 Set right order
`utilities/fasta2paml_get_order.py`, see --help for args. Folder with .order files can be empty, files with right order will be recorded to that folder. 
#### 4.3 Preprocessing, convert fasta to paml format  
##### 4.3.1 fasta2paml.py
Skip if right order was set in the step 4.2.
The script (fasta2paml.py) consists of two stages:  
1. Converting fasta format nucleotide codon sequences (from input directory) to philip-sequential format (to output 
directory)  
2. Converting philip-sequential format to specific philip format required by PAML:  
In resulting out_dir:  directory of name "group_id" with folders "file_name" with file_name.phy file for PAML.
For example: 
```bash
$ cd /abspath/tothe/nuc_out_prank/  
$ ls */  
12/:            23/:           12345/:
1.fasta         4055.fasta     2031.fasta   
                3010.fasta     2.fasta
```
Result:
```bash
$ cd out_dir  
$ ls */  
12/:            23/:            12345/:
1/:             4055/:          2031/:   
1.phy           4055.phy        2031.phy
                3010/:          2/:
                3010.phy        2.phy    
```
In --i and out --o folders can be the same. Making backups is necessary for further analysis.  
`python pipeline/fasta2paml.py --i /abspath/tothe/nuc_out_prank/ --o /abspath/tothe/nuc_out_prank/  
--species 8 --group 6`  
See 'fasta2paml.log' in working directory.  
###### 4.3.1.1 Test right order manually
Test PAML (codeml) for knowing right order for sequences to exclude PAML's errors (Some reference to tree file format from [PAML manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf): "The species can be represented using either their names or their indexes corresponding to the order of their occurrences in the sequence data file.", but there may be some nuances here).
Test can be performed with launching the one ratio PAML model with script 'paml_one_ratio_model.py',  
see --help for arguments: option --e can be skipped if use codeml from biopython (Bio.Phylo.PAML), use option '--rework y' if want to overwrite existing paml files. This script writes .ctl file and launches one ratio model.<br />
`python pipeline/paml_one_ratio_model.py --i /abspath/tothe/nuc_out_prank/ --tree /abspath/folder_trees/ --threads 22`    
See "paml_one_ratio.log", further testing may be continued in separate item's folders just from command line.
##### 4.3.2 Ordering
- After having known about the right orders, .order files should be placed to separate folder. Name of .order file should be the same as name of species folder ('12' -> '12.order');<br />
- Launch fasta2paml_ordering.py (can be launched on the backup folder), right order will be set automaticaly:<br />
`python pipeline/fasta2paml_ordering.py --i /abspath/tothe/nuc_out_prank/ --order /abspath/folder_orders/`<br />`--o /abspath/tothe/nuc_out_prank/ 
--species 8 --group 6` 
See 'fasta2paml_ordering.log' in working directory.  
#### 4.6 One ratio model  
See launch example above in the step 4.3.1.1.
#### 4.7 SWAMP masking
Sliding window approach SWAMP to mask regions of the alignment with excessive amino acid changes.
- Construct branchcodes for every species group:  
  - required tree view for automatic build branchcode (items are separated by spaces): '(1, ((2, 3), ((6), 8, 4)));'
  - if you have folder with marked trees for paml, you can clean it from label (#1) and insert spaces with sed stream editor:<br />
  `sed -i 's/ #1//' *` <br />
  `sed -i 's/,/, /g' *`
  - `python pipeline/construct_branchcodes.py --i /abspath/tothe/nuc_out_prank/ --t /abspath/folder_trees_clean/`<br />`--b /abspath/folder_for_branchcodes/`
- Launch SWAMP:
  - "swamp_script.py" for python3 environment, swamp_script_py2.py for python2 envoronment;
  -  use modified version of SWAMP executable GScour/"SWAMP_ordered.py" to conserve right order.<br />
`python pipeline/swamp_script_py2.py -e /GScour/SWAMP_ordered.py -i /abspath/tothe/nuc_out_prank/`<br />` -b /abspath/tothe/branchcodes/ -t 2 -w 20`<br />
See stdout and 'swamp_log.log'. <br />
Use global variable 'target_dict' in swamp_script.py if need to run on individual files:<br />
`target_dict[species_folder] = [item_folder1, item_folder2...]`
#### 4.8 Perform maximum likelihood (ML) dN/dS analysis to infer positive selection of genes and codons, using codeml from the PAML software package.
__Branch-site model__  
There are two hypothesis:  
```
H0 (The null model for Branch-site model A):  
    Model A1: model = 2, NSsites = 2, fix_omega = 1, omega = 1  
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated  
    kappa = 2   * initial or fixed kappa  
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate  
    omega = 1   * initial or fixed omega, for codons or codon-based AAs  
H1 (Alternative model, Model A: model = 2, NSsites = 2, fix_omega = 0 ):  
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated  
    kappa = 2   * initial or fixed kappa  
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate  
    omega = 1   * initial or fixed omega, for codons or codon-based AAs  
```    
From readme of paml example lysozymeLarge.ctl:
Alternative hypothesis (branch site model A, with w2 estimated):  
```
model = 2    NSsites = 2   fix_omega = 0   omega = 1.5 (__or any value > 1__)
```
> As the branch-site model is known to cause computational difficulties for the numer-  
ical optimization algorithm in PAML, each analysis is conducted three times with  
different starting values to ensure that the global peak is found (*Statistical Properties of the Branch-Site Test of Positive
Selection, Ziheng Yang and Mario dos Reis*)

- PAML launch
For masking files (after SWAMP) launch, for example:  
`python pipeline/masked_paml_branch_site_model.py --timeout 1000 --i /abspath/tothe/for_paml/ --tree /abspath/folder_trees/ --threads 64 --rework y` 

For files without masking:  
`python pipeline/paml_branch_site_model.py --timeout 1000 --i /abspath/tothe/for_paml/ --tree /abspath/folder_trees/ --threads 64 --rework y` 

See help for args.

- Additional check<br />
Run `utilities/count_correct_rst_files.py` to check number of correct auxiliary files required for paml analysis.<br />

### 5. Analysing PAML's results  
#### 5.1 Tabulation
`python pipeline/paml_out_analysis.py (paml_out_analysis_masked.py) --i /abspath/tothe/for_paml/` <br />`--log /abspath/tothe/nuc_out_folder/ --required 6`<br />
Results will be written to /abspath/tothe/for_paml/common_sheet.xlsx, also will be written in every species folder to file named 'name_of_species_folder.result' and summary to stdout.
* common_sheet.xlsx, sheet for every species group:<br />
dN/dS (w) for site classes (K=4) (see [PAML manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf))<br />
site class 0 1 2a 2b <br />
proportion (proportion of sites that have those omega values, for each site class) <br />
background w (omega values on the background, for each site class) <br />
foreground w (omega values on the foreground, for each site class) <br />
P-value, likelihood from positive vs lokelihood from neutral model. <br />
* summary sheet with genes under positive selection (p-value < 0.05). P-value can be adjusted here: <br />
https://github.com/ctlab/GScour/blob/a2145a12a754e94b6306dce69bfcd7b173d8a898/pipeline/paml_out_analysis.py#L131  
#### 5.2 Bonferroni and FDR correction
* Take into account the information from (Yang, 2007):
"The branch-site test requires a priori specification of the foreground branches. When multiple branches on the tree are tested for positive selection using the same data set, a correction for multiple testing is required (Anisimova and Yang 2007). A simple and slightly conservative procedure is Bonferroni's correction, which means that the individual test for any branch is considered significant at the level α only if the p-value is <α/m, where m is the number of branches being tested using the same data."
* Use "adjust_pvalue.py" script, see --help for args and adjust the denominator if necessarily in lines:
https://github.com/ctlab/GScour/blob/ea7787697646523899088ef94ccb1e5f3cb0e935/utilities/adjust_pvalue.py#L24
https://github.com/ctlab/GScour/blob/ea7787697646523899088ef94ccb1e5f3cb0e935/utilities/adjust_pvalue.py#L26
#### 5.3 Assembling results
Script "assembling_results.py" collects final data from sheets "summary" for every result common sheet: concatenate and remove duplicates. 

