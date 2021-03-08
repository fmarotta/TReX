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
each gene. By default, only promoters are used. We recommend to split
the files by chromosome and work on one chromosome at a time.

**Input**

    1. A GTF or GFF3 file with all the transcripts annotation

    2. (Optional) A BED file with the enhancers

**Command**

```
    Rscript make_regreg.R \
        transcripts.gtf \           # The GTF or GFF3 file with the annotations
        --upstream=1500 \           # Promoter start upstream of the TSS
        --downstream=500 \          # Promoter end downstream of the TSS
        --chromosome=chr7 \         # Consider only the specified chromosome
        --threads=1 \               # Number of cores for parallel computations
        --output=chr7.regreg.bed    # Path to the output file
```

If one wishes to include enhancers in addition to promoters, they
can separately create or download a BED file with all the enhancers
coordinates and supply this file to `make_regreg.R` with the option
`--enhancers=<PATH_TO_FILE>`. Then, to each gene it will be assigned the
closest enhancer.

The output of `make_regreg.R` will be a BED file with the coordinates
and the names of all the regulatory regions associated to each gene.

### 2. Computing the total binding affinity

Having defined the regulatory regions, it is now time to compute the
total binding affinities of transcription factors for the regions.

**Input**

1. BED file with the regulatory regions (obtained in the previous step)

2. VCF file with the genetic variants of each sample

3. Position weight matrices of all transcription factors in TRANSFAC
format

4. FASTA file with the reference sequence of the genome

5. (Optional) Background letter frequencies of the genome

**Remarks**

We recommend some preprocessing of the VCF before using it with TReX.
For example, it may be useful to remove samples or variants with too
many missing values, or to keep only variants with a MAF of at least
0.01. `plink2` is an excellent resource for this type of processing.

The 'Background letter frequencies' file should be a simple file
with one line and four columns: each column should contain the frequency
of A, C, G, and T, respectively, in the reference genome.

**Command**

```
bash compute_tba.sh \
    --bed chr7.regreg.bed \         # The regulatory regions
    --vcf chr7.vcf.gz \             # The VCF
    --pwm human_motifs.transfac \   # The position weight matrices
    --ref chr7.fa.gz \              # The reference sequence
    --bkg bgk_frequencies.txt \     # The background frequencies
    --threads 1 \                   # Number of threads for parallel computations
    --outdir tba/chr7/              # Output directory
```

The computation should take around 1-20 hours for each chromosome
(parallelisation won't help much). The two main outputs that will be
created by `compute_tba.sh` are *tba.tsv* and *vcfrider_snps.tsv*,
the former containing the total binding affinity values for each
gene-sample-transcriptionfactor triplet.

### 3. Fitting the expression-predicting model

Now, using a data set where both genotypes and gene expression are
available (*e.g.* GTEx, ROS-MAP, ...), we train the regression models
that predict gene expression from the total binding affinities.

**Input**

1. Total binding affinities file (obtained in the previous step)

1. Gene expression file

1. (Optional) List of samples ID

**Remarks**

We recommend some preprocessing for the gene expression
file; best practices are described at the [GTEx
portal](https://www.gtexportal.org/home/documentationPage#staticTextAnal
ysisMethods). In particular, we recommend normalization and
inverse-normal transformation of the expression values. As TReX does
not (yet) support custom covariates in the regression model, we also
recommend to not work with the expression directly, but rather with the
residuals of the model `EXPRESSION ~ COVARIATES` (The covariates can be
downloaded from the GTEx portal).

**Command**

```
Rscript fit_model.R \
    expr_residuals.tsv \                # Path to the expression file
    tba/chr7/tba.tsv \                  # Path to the TBA file
    --folds 1 \                         # Number of folds for performance evaluation
    --nested-folds 10 \                 # Number of folds for parameter tuning
    --min-R2 0.1 \                      # Minimum performance R^2 to consider
    --train-samples samples_ID.txt \    # List of samples
    --threads 1 \                       # Number of threads for parallel computations
    --outdir fit/chr7/                  # Output directory
```

**Tips**

Since we want to evaluate the performance of a model and tune some
parameters at the same time, we must use a 'nested cross-validation.'
The `folds` argument specifies the number of folds in the outer loop; in
the existing literature, this is set to 5. If you don't want to evaluate
the performance (it is time consuming after all), pass 1 to this
argument. The `nested-folds` argument specifies the number of folds in the inner
loop; by default, it is 10. After the cross-validation, the model with
the best parameters is fit on the whole data. As a measure of the
performance of the model, we use the R^2^ (the square between the true
and the predicted expression values).

When you choose to evaluate the performance (thereby setting the `folds`
argument to something greater than 2), you can decide to not bother
fitting the model if the performance was too low. The `min-R2` argument
lets you choose the minimum cross-validation R^2^ after which the model
will be fitted. Recall that there are two steps: first the
cross-validation evaluates the performance and finds the optimal
parameters, then the model is fitted on the whole data. If, after the
cross-validation, the performance was low (*i.e.* the R^2^ was low), it
means that predicting the expression of this gene is too difficult for
the model, so it may be best to discard that gene entirely.

The *samples_ID.txt* file should be a simple one-column file with one
sample ID for each line. The individuals listed in this file will be
used for the training. You may want to use this argument to save some
individuals for the validation or the testing of the model. By default,
all individuals are used.

This is likely to be the slowest step of the pipeline, but it can be
parallelised at will.

You don't have time to train these models? E-mail the authors of the
paper: they will be more than happy to provide you with the pre-trained
models. Maybe they will even be uploaded to a public server.

### 4. Predicting gene expression in new samples

The models were trained in a data set where both genotypes and gene
expressions were available. Now we apply the model to a data set where
genotypes and phenotypes are available. But first, we have to predict
the gene expression in those individuals whose genotypes are available.
The first step here is to compute the total binding affinities in those
samples. The procedure is the same as in section [2. Computing the total
binding affinity](#2.-Computing-the-total-binding-affinity), except that
now the VCF file refers to the new data set.

