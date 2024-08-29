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
    -a, --adaptors: The lengths of your 5' and 3' adaptors separated by a comma (e.g., 5,3).
    -t, --title: Title of the plot. Default is no title.
    -d, --difference: Calculate the difference between two input files (second input - first input).
    -b, --bidirectional: Produce a difference plot with directionality maintained (color bar normalized from [-0.5, 0.5]).
    -r, --resolution: Resolution of the final image in DPI. Default is 300.

### Example Usage

Regular heatmap

    python tRNAheatmap.py -i WT-IVT_c0.vcf -o yeast_wildtype_heatmap.png -l tRNA_lengths.txt -a 18,6 -t "Wildtype"

Difference heatmap

    python tRNAheatmap.py -i WT-IVT_c0.vcf pus4-IVT_c0.vcf -o yeast_pus4d-WT_heatmap.png -l tRNA_lengths.txt -a 18,6 -t "pus4d - WT" -d


## Requirements
- Python 3.x
- NumPy
- Matplotlib
