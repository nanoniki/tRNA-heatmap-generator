# tRNA Posterior Probability Heatmap Generator
## Overview
This Python script generates a heatmap of posterior probabilities along tRNA isoacceptors sequenced with Nanopore technology. It takes input files generated with marginCaller and produces a graphical representation of the posterior probabilities for each tRNA isoacceptor.

![alt text](https://github.com/nanoniki/tRNA-heatmap-generator/blob/main/Examples/WT_structure_example.png)

## Features
- Supports single or paired input files for comparison.
- Calculates the difference between two input files if specified.
- Provides options for customizing the plot title, resolution, and adaptors.
- Handles bidirectional difference plots to maintain directionality.
- Generates high-quality heatmaps in various output formats (e.g., PDF, PNG, JPEG).

## Usage
### Input Files
The input files should be generated with marginCaller using the --threshold=0 setting to ensure that every position is in the output. Though not required, we recommend using marginAlign EM training with synthetic or IVT data as your error model for the most accurate results. If you're comparing two files to calculate the difference, ensure that they have the same format and positional information.

### Command-Line Arguments
    -i, --input: One or two input pileup .vcf files generated with marginCaller. If --difference is used, provide two files.
    -o, --output: Path to the output location with the desired file extension (e.g., heatmap.png).
    -l, --lengths: Path to a tab-separated file containing the name of each tRNA and their respective lengths with adaptors. The naming conventions used for tRNA in this file will appear as y-axis labels.
    -s, --structure: Path to a tab separated file containing the name of each tRNA and the respective coordinates of their variable positions.
    -a, --adaptors: The lengths of your 5' and 3' adaptors separated by a comma (e.g., 5,3).
    -t, --title: Title of the plot. Default is no title.
    -d, --difference: Calculate the difference between two input files (second input - first input).
    -b, --bidirectional: Produce a difference plot with directionality maintained (color bar normalized from [-0.5, 0.5]).
    -r, --resolution: Resolution of the final image in DPI. Default is 300.

### Example Usage

Regular heatmap

    python tRNAheatmap.py -i WT-IVT_c0.vcf -o yeast_wildtype_heatmap.png -l tRNA_lengths.txt -s variable.txt -a 18,6 -t "Wildtype"

Difference heatmap

    python tRNAheatmap.py -i WT-IVT_c0.vcf pus4-IVT_c0.vcf -o yeast_pus4d-WT_heatmap.png -l tRNA_lengths.txt -s variable.txt -a 18,6 -t "pus4d - WT" -d


## Requirements
- Python 3.x
- NumPy
- Matplotlib

## Publication Data
### Authors Table 1
Sample  refers   to   nomenclature   used   in   Shaw   and   Thomas  et al.  2024   (PMID:
39340295). Fast5 File  and  Fastq File  refer to file names that are publicly accessible at the
European Nucleotide Archive https://www.ebi.ac.uk/ena. The study accession number
is PRJEB68201.

| Sample                       | Fast5 File                                    | Fastq File                                    |
|------------------------------|-----------------------------------------------|-----------------------------------------------|
| Wild-type Replicate 1         | WT_tRNAs_rep1_7_22_21_fast5s.tar.gz           | WT_tRNAs_rep1_7_22_noU.fastq.gz               |
| Wild-type Replicate 2         | WT_tRNAs_rep2_8_6_21_fast5s.tar.gz            | WT_tRNAs_rep2_8_6_21_noU.fastq.gz             |
| Wild-type Replicate 3         | WT_tRNAs_rep3_9_17_21_fast5s.tar.gz           | WT_tRNAs_rep3_9_17_21_noU.fastq.gz            |
| Wild-type Replicate 4         | WT_tRNAs_rep4_4_4_22_fast5s.tar.gz            | WT_tRNAs_rep4_4_4_22_noU.fastq.gz             |
| pus4Δ Replicate 1             | pus4_tRNAs_rep1_06_30_21_fast5s.tar.gz        | pus4_tRNAs_rep1_06_30_21_noU.fastq.gz         |
| pus4Δ Replicate 2             | pus4_tRNAs_rep2_07_19_21_fast5s.tar.gz        | pus4_tRNAs_rep2_07_19_21_noU.fastq.gz         |
| pus4Δ Replicate 3             | pus4_tRNAs_rep3_10_26_21_fast5s.tar.gz        | pus4_tRNAs_rep3_10_26_21_noU.fastq.gz         |
| trm6Δ Replicate 1             | gcd10_tRNAs_rep1_1_30_23_fast5s.tar.gz        | gcd10_tRNAs_rep1_1_30_23_noU.fastq.gz         |
| trm6Δ Replicate 2             | gcd10_tRNAs_rep2_2_28_23_fast5s.tar.gz        | gcd10_tRNAs_rep2_2_28_23_noU.fastq.gz         |
| trm6Δ Replicate 3             | gcd10_tRNAs_rep3_2_28_23_fast5s.tar.gz        | gcd10_tRNAs_rep3_2_28_23_noU.fastq.gz         |
| trm2Δ Replicate 1             | trm2_tRNAs_rep1_12_20_22_fast5s.tar.gz        | trm2_tRNAs_rep1_12_20_22_noU.fastq.gz         |
| trm2Δ Replicate 2             | trm2_tRNAs_rep2_6_30_22_fast5s.tar.gz         | trm2_tRNAs_rep2_6_30_22_noU.fastq.gz          |
| trm2Δ Replicate 3             | trm2_tRNAs_rep3_6_24_22_fast5s.tar.gz         | trm2_tRNAs_rep3_6_24_22_noU.fastq.gz          |
| pus4Δtrm2Δ Replicate 1        | pus4_trm2_tRNAs_rep1_12_20_22_fast5s.tar.gz   | pus4_trm2_tRNAs_rep1_12_20_22_noU.fastq.gz    |
| pus4Δtrm2Δ Replicate 2        | pus4_trm2_tRNAs_rep2_8_11_22_fast5s.tar.gz    | pus4_trm2_tRNAs_rep2_8_11_22_noU.fastq.gz     |
| pus4Δtrm2Δ Replicate 3        | pus4_trm2_tRNAs_rep3_12_20_22_fast5s.tar.gz   | pus4_trm2_tRNAs_rep3_12_20_22_noU.fastq.gz    |
| PUS4 R286K Replicate 1        | pus4_catdead_tRNAs_rep1_10_5_21_fast5s.tar.gz | pus4_catdead_tRNAs_rep1_10_5_21_noU.fastq.gz  |
| PUS4 R286K Replicate 2        | pus4_catdead_tRNAs_rep2_2_9_22_fast5s.tar.gz  | pus4_catdead_tRNAs_rep2_2_9_22_noU.fastq.gz   |
| PUS4 R286K Replicate 3        | pus4_catdead_tRNAs_rep3_1_16_22_fast5s.tar.gz | pus4_catdead_tRNAs_rep3_1_16_22_noU.fastq.gz  |
| IVT                          | 07_19_22_R941_yeastIVTtRNA_fast5s.tar.gz      | IVT_08-16-22_NoU.fastq.gz                     |
| Wild-type Saturated Replicate 1| WT_saturated_tRNAs_rep1_11_8_21_fast5s.tar.gz | WT_saturated_tRNAs_rep1_11_8_21_noU.fastq.gz  |
| Wild-type Saturated Replicate 2| WT_saturated_tRNAs_rep2_11_8_21_fast5s.tar.gz | WT_saturated_tRNAs_rep2_11_8_21_noU.fastq.gz  |
| Wild-type Saturated Replicate 3| WT_saturated_tRNAs_rep3_2_1_22_fast5s.tar.gz  | WT_saturated_tRNAs_rep3_2_1_22_noU.fastq.gz   |
| trm6Δ Saturated               | gcd10_saturated_tRNAs_rep1_fast5s.tar.gz      | gcd10_saturated_tRNAs_rep1_noU.fastq.gz       |

