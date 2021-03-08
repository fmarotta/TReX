# TReX

A suite of scripts to predict gene expression from the DNA sequence
using position weigth matrices, and perform a transcriptome-wide
association study (TWAS).

# Installation

1. Clone the repository (`git clone git@github.com:fmarotta/TReX.git`)

2. Install `vcf_rider` (follow the instructions at the [project's
repository](https://github.com/vodkatad/vcf_rider))

3. Install `blort` (follow the instructions at the [project's
repository](https://github.com/fmarotta/blort))

4. Launch R and install the required packages from either CRAN or
Bioconductor, as appropriate:
    * GenomicFeatures
    * docopt
    * fplyr

# Usage

Transcriptome-wide association studies can be performed in two modes:
genotype-level, when all the genotypes of the individuals are available;
and summary-level, when only summary statistics are available. While
TReX supports both of these modes, the two pipelines are somewhat
different and will be described separately.

## Genotype-level TWAS

### 1. Define the regulatory regions

TReX uses the total binding affinities of transcription factors for the
regulatory regions of a gene in order to predict the gene's expression;
The first step of the analysis is to define the regulatory regions for
each gene. By default, only promoters are used.

*Input*
    1. A GTF or GFF3 file with all the transcripts annotation

*Command*
```
    Rscript make_regreg.R \
        --upstream=<N> \        # Promoter start upstream of the TSS [default: 1500]
        --downstream=<N> \      # Promoter end downstream of the TSS [default: 500]
        --chromosome=<chrN> \   # Consider only the specified chromosome [default: all]
        --output=<FILE> \       # Path to the output file [default: STDOUT]
        --threads=<N>           # Number of cores for parallel computations [default: 1]
```

### 2. Computing the total binding affinity

The first step is to compute the total binding affinity of the
regulatory regiosn
