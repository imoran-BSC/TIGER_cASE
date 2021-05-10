# TIGER cASE analysis repository
Combined ASE expression scripts used in the Alonso _et al._ TIGER paper.

The results generated by the scripts available in this repository, plus 
all other results generated in this study, can be explored in the
[TIGER Data Portal](http://tiger.bsc.es), alongside with other
relevant genomic, epigenetic and transcriptomic data.

## Preface
This repository contains the C++ source code for four scripts that
were used in different analyses in the study.

Each script can be compiled to a binary executable file by using

```g++ -Wall -o stepN_binary step_N.cpp ../utils/allelic_imba_common.cpp```

Each binary expects a ```config.ini``` file residing in the same directory,
containing the necessary arguments for the executions, in their expected 
order. Example config files are already included in this repository, with 
a brief description of  the type of the expected arguments, and of 
their particular formatting.
Prior to executing the binary, please _make sure that all the necessary 
arguments in ```config.ini``` are correct_: that the order of the arguments
was not altered, and that they refer to files that exist and that
are formatted as expected. Otherwise, the execution will terminate unexpectedly.

The resulting binary can then be executed in a unix terminal as follows:

```./stepN_binary```

These scripts are numbered to indicate their execution order in the 
analysis workflow. However, these don't comprise a full pipeline, 
since some can be executed in a standalone manner, or some require
pre- or post-processing of the data.

The code is provided as is, for informational purposes, and without any 
warranties. It has only been tested in a limited number of environments.

## 01 Artificial Reads
This script reads a list of SNPs, a gtf gene annotation and a fasta genome, 
and generates all possible reads overlapping the given SNPs, including spliced
reads across splice junctions. The resulting set of artificial reads includes 
the same number of  reads containing the reference and alternate alleles 
for all SNPs. It also outputs
a histogram of the SNP density in each read position, and a list of regions
discarded due to excessive density of polymorphisms.

A set of artificial RNA-seq reads generated with this code was used to 
test for alignment allelic biases in all possible reporter variants. Given 
that the same number of reads are generated containing the reference and 
alternate alleles, any deviation form the expected 50% allelic ratio can be
attributed to problems with the alignment strategy, and discarded
from downstream analyses.

## 02 Read Merger
Reads the indexed and sorted .bam files outputed from the ENHANCED 
and MASKED alignments, and outputs a MERGED standard and non-clonal
.bam files.

This script _requires ```samtools``` to be installed and system-wide
available_, and it has been tested with version 1.7. Other versions have 
not been tested, and if they have different syntax requirements 
execution will unexpectedly fail.
The system calls use unix-specific syntax, so this code is not expected
to work in other environments.

## 03 Filter and Unify Alignments
Given a list of sample names as input, this script reads all allelic 
imbalance quantification.txt files, obtained using ```samtools mpileup``` and 
```ComputePileupFreqs.pl``` on the output files of step 02,
as explained in the Online Methods section. 
It also requires a .vcf genotype file, with the genotypes 
of all samples of interest in columns, in any order. This script then
filters all input data and creates a unified output file with all the 
allelic imbalance information, along with other files on mean allelic biases,
or Benjamini Hochberg significance thresholds.

The resulting outputs are intermediary files necessary for downstream 
analyses.

## 04 Permuted Reporter Imbalance
This script reads an Allelic_Imba_trimmed.txt file, outtputed in step 02,
calculates the empiric Z-score distribution, then randomises
the Het variant read counts N times, and creates an expected null Z-score 
distribution. The counts are randomised within 5 gene expression bins, 
one for SNPs with a median read coverage of zero, and four more in increasing 
quartiles for the remainder of the data, to control for possible biases
induced by gene expression level or read coverage depth.

The Z-score distribution obtained from this script is used 
in the downstream analyses as the null distribution, _i.e._ the expected
distribution of allelic imbalance Z-scores that would be obtained if 
all allelic biases were the product of biological stochasticity.

## 05 Find Candidate Variants
This script reads an Allelic_Imba_trimmed.txt file outputted in step 02,
the control Z-score distribution outputted in step 03,
and the full genotype of the given chr. It also reequires PLINK 
preprocessing of the sample genotypes, to obtain linkage disequilibrium, 
allele frequency and Hardy-Weinberg information, split by chromosome.

For each reporter SNP, this script 
calculates the reporter Z-score (outputted in Allelic_Imba_reporters_chrN.txt).
A larger deviation from a zero value indicates a stronger allelic bias in
the expression of the reporter, with positive values indicating a preference
for the reference allele, and negative for the alternate.

Additionally, this script also calculates the Z-scores for each candidate 
variant around each reporter, within a TAD+leeway (outputted in Allelic_Imba_CSet_chrN.txt). 
It separates the samples in two populations, those Het and those Hom for 
the candidate, as described in the Online Methods section. 
If the Het Z-score is stronger than the reporter,
and the Hom Z-score is non significant, the variant in question is 
identified as a putative candidate for the _cis_-regulatory effect.


