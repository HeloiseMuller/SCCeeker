# SCCeeker

## Introduction
SCCeeker has been developed to annote SCCmec and SCCmec-like elements.

## Citation
Please cite our paper "Diversification of the staphylococcal cassette chromosome through distinct mechanisms of horizontal transfer"  by Héloïse Muller and Cheryl P. Andam (Publication to come)

## Requirements
- [R](https://cran.r-project.org) 4.3+ with the following packages:
  - [data.table](https://cran.r-project.org/web/packages/data.table/) 1.16.2
  - [stringi](https://cran.r-project.org/web/packages/stringi/) 1.8.4
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/) 1.1.4
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/) 3.5.1
  - [stringr] 1.5.1
  - [gridExtra] 2.3
  - [ggpubr] 0.6.0

  (If not found, these packages are installed automatically by the pipeline.) 

- [ncbi blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.14.0
- [python3]

The pipeline was not tested with other versions of the above programs, but other versions probably work.  

## Installation
In a bash-compatible terminal that can execute git, paste
```
 git clone https://github.com/HeloiseMuller/SCCeeker.git
 cd SCCeeker/
 unzip DB_20171117.zip 
```

## Get files ready to run SCCeeker
SCCeeker needs a two columns table containing the absolute path of each genome you want to scan for SCCmec elements and the pattern to use in all outputs for each file. See ` test/example_input.tbl ` as an example. 

## Test SCCeeker
1) Replace the PATH found in ` test/example_input.tbl `  by the correct path
2) Run the next two command lines,  by replacing paths:
```
bash SCCeeker_severalGenomes.sh -d $ABSOLUTE_PATH_SCCeeker/SCCeeker/DB_20171117/  -f $ABSOLUTE_PATH_SCCeeker/test/example_input.tbl -k extended -o $ANY_PATH 
Rscript SCCeeker.R -f $ANY_PATH /cat_all5_blast.out -w whole_cassette_SCCmec_database_EXTENDED_20171117 -k $ANY_PATH /cat_all5_kmer.tsv 
```
where `$ABSOLUTE_PATH_SCCeeker` is the path where you cloned SCCeeker, $ANY_PATH  is the path where you want the ouputs to be saved at.  

To get all options, runs:
```
bash SCCeeker_severalGenomes.sh -h
Rscript SCCeeker.R -h
```

## Output description of SCCeeker.sh
 
`cat_allX_blast.out` where X is the total number of analyzed assemblies is the blast output.  
`cat_allX_kmer.tsv` where X is the total number of analyzed assemblies is the output of the k-mer approach if the option was activated.

## Output description of SCCeeker.R

`SCCeeker_summary_perRegion.out` summarizes all candidate regions identified in the dataset. One assembly can have several candidate regions.  
`SCCeeker_summary_perFile.out` is a summary for each assembly. Are reported the presence of a mec gene, the number of candidate regions, the type of cassette if any region could be typed, and the same but reported the type only if it was validated by the k-mer approach.   
Many intermediate filed can be found in the directory detailed_outputs/  
Some figures summarizing the results can be found in the directory Figures/  

## Content of SCCeeker

`DB_20171117/` was downloaded directory from SCCmecFinder. It contains :
 - single_genes_database_20171117.fasta: genes of interest to identify and type cassettes
 - mec_database_20171117.fasta: fasta to differentiate betwee mec class C1 or mec class C2
 - whole_cassette_SCCmec_database_REFERENCE_20171117.fasta: whole length cassettes
 - whole_cassette_SCCmec_database_EXTENDED_20171117.fasta: extended version of the whole length cassettes
 - selectDB: file that described which of the above database to use. To modify only of one wants to update the database. By default, SCCeeker will use the EXTENDED database.
 - template_db: database used for the k-mer approach 

`test/` contains inputs to test SCCeeker

`SCCeeker.sh` runs the first step of SCCeeker (blast and an optional k-mer approach) on one genome.  
Alternatively, `SCCeeker_severalGenomes.sh` runs the first step of SCCeeker on several genomes.  
`findtemplate.py` is the k-mer approach ran by the first step of SCCeeker. This script comes from SCCmecFinder.  
`SCCeeker.R` is the second step of SCCeeker that identify candidate regions and attempt to type them.   
`typing.R` contains the definition of ccr complexes, mec complexes and the 15 approved typed. It is read by `SCCeeker.R`. One might want to update it when more typed will be approved. Otherwise, no need to touch it.  