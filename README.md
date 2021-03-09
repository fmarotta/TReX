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

## Model training

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
Consider also the imputation of missing genotypes, especially if the
genotypes don't come from whole-genome sequencing.

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
portal](https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods).
In particular, we recommend normalization and
inverse-normal transformation of the expression values. As TReX does not
support custom covariates in the regression model, we also recommend to
not work with the expression directly, but rather with the residuals of
the model `EXPRESSION ~ COVARIATES` (The covariates can be downloaded
from the GTEx portal).

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
    --outdir training/chr7/             # Output directory
```

**Tips**

Since we want to evaluate the performance of a model and tune some
parameters at the same time, we must use a 'nested cross-validation.'
The `folds` argument specifies the number of folds in the outer loop;
in the existing literature, this is set to 5. If you don't want to
evaluate the performance (it is time consuming after all), pass 1 to
this argument. The `nested-folds` argument specifies the number of folds
in the inner loop; by default, it is 10. After the cross-validation, the
model with the best parameters is fit on the whole data. As a measure
of the performance of the model, we use the R<sup>2</sup> (the square
between the true and the predicted expression values).

When you choose to evaluate the performance (thereby setting the `folds`
argument to something greater than 2), you can decide to not bother
fitting the model if the performance was too low. The `min-R2` argument
lets you choose the minimum cross-validation R<sup>2</sup> after which
the model will be fitted. Recall that there are two steps: first the
cross-validation evaluates the performance and finds the optimal
parameters, then the model is fitted on the whole data. If, after the
cross-validation, the performance was low (*i.e.* the R<sup>2</sup>
was low), it means that predicting the expression of this gene is
too difficult for the model, so it may be best to discard that gene
entirely.

The *samples_ID.txt* file should be a simple one-column file with one
sample ID for each line. The individuals listed in this file will be
used for the training. You may want to use this argument to save some
individuals for the validation or the testing of the model. By default,
all individuals are used.

This is likely to be the slowest step of the pipeline, but it can be
parallelised at will. Nevertheless, if you cannot train these models,
email the authors of the paper: they will be more than happy to provide
the pre-trained models. Maybe the models will be even uploaded to a
public server.

## Genotype-level TWAS

### 1. Predicting gene expression in new samples

After fitting the model (or downloading the precomputed weights), it
is time to apply it to a data set where genotypes and phenotypes (but
not expressions) are available. But first, we have to predict the gene
expression in those individuals whose genotypes are available. The first
step here is to compute the total binding affinities in those samples.
The procedure is the same as in section [2. Computing the total binding
affinity](#2-computing-the-total-binding-affinity), except that now the
VCF file refers to the new data set.

**Input**

1. BED file with the regulatory regions (obtained at the beginning)

2. VCF file with the genetic variants of each sample in the new data
set

3. Position weight matrices of all transcription factors in TRANSFAC
format

4. FASTA file with the reference sequence of the genome

5. (Optional) Background letter frequencies of the genome

**Remarks**

If you have just downloaded the precomputed weights, you
should still define the regulatory regions as described
[above](#1-define-the-regulatory-regions).

We recommend to apply the same preprocessing to the VCF in the training
data set and in the prediction data set. Consider also the imputation of
missing genotypes.

Make sure that the genome assembly versions are the same between the two
data set, otherwise you will have to convert them. For instance, both
data set should have hg38 coordinates.

**Command**

```
bash compute_tba.sh \
    --bed chr7.regreg.bed \             # The regulatory regions
    --vcf chr7.phenodataset.vcf.gz \    # The VCF
    --pwm human_motifs.transfac \       # The position weight matrices
    --ref chr7.fa.gz \                  # The reference sequence
    --bkg bgk_frequencies.txt \         # The background frequencies
    --threads 1 \                       # Number of threads for parallel computations
    --outdir tba/genopheno/chr7/        # Output directory
```

Now that we have the predictors, we can use the training weights to
predict gene expression.

**Input**

1. Total binding affinities file (obtained in the previous step)

1. R object with the model fit (obtained in the training phase)

1. (Optional) List of samples ID

**Command**

```
Rscript ../R/predict_trex.R \
    training/chr7/fit/trex_fit.Rds \            # The fit object
    tba/genopheno/chr7/tba.tsv \                # The TBA file
    --test-samples samples/genopheno/list.txt \ # The list of samples
    --threads 1 \                               # Number of cores for parallel computations
    --outdir prediction/chr7/                   # Output directory
```

**Tip**

If the true expression for this data set is available, some statistics
about the accuracy of the prediction can be generated by passing the
true expression file path with the `--expression` argument.

### 2. Performing the genotype-level TWAS

Finally, we associate the predicted expression to the phenotype. The
association test is made with either linear or logistic regression,
according to whether the phenotype is quantitative or binary. The output
of this step will be a file containing, for each gene, an estimate of
the association (the coefficient of the regression), its standard error,
its Z-score, and its p-value.

**Input**

1. Predicted expression file (generated in the previous step)

2. Phenotype file

3. Covariates file

**Remarks**

The phenotype file can actually contain multiple phenotypes. It should
contain the sample IDs in the first column and all the phenotypes in
the subsequent columns (one sample per line). TReX will perform one
association test for each phenotype, after merging the covariates.

The covariates file should have the sample IDs in the first column
and all the covariates (*e.g.* age, sex, principal components of the
genoype...) in the next columns.

**Command**

```
Rscript genotype_twas.R \
    phenotypes/alzheimer.tsv \                  # The phenotypes file
    prediction/chr7/trex_pred.tsv \             # The predicted expression
    --covar phenotypes/alzheimer_covar.tsv \    # The covariates
    --threads 1 \                               # Number of cores for parallel computations
    --outdir twas/individual/alzheimer/chr7/    # Output directory
```

**Tip**

Use the `--pheno-quantile-normalize` to apply an inverse normal
transformation to the quantitative phenotypes before performing the
regression.

## Summary-level TWAS

### 1. Compute the delta-TBA

In this step the effect of each SNP on the total binding affinity will
be computed. Later, the change in TBA caused by each SNP will be
combined with the change in gene expresson caused by the TBA, to
ultimately find the change in gene expression caused by each SNP. Then,
the change in gene expression caused by each SNP (through the TBA) will
be combined with the GWAS Z-score of that SNP to find associations
between SNPs and phenotypes.

This step must be performed on a reference data set where the genotypes
are available. We recommend to use, if possible, the same data set used
for the training. If you have downloaded the precomputed weights and did
not perform the training, ask the authors for the precomputed delta-TBA
as well.

**Input**

1. BED file with the regulatory regions (see
[above](#1-define-the-regulatory-regions))

2. VCF file with the genetic variants of each sample in the training
data set

3. Total binding affinities file (see
[above](#2-computing-the-total-binding-affinity))

3. vcf\_rider's SNPs file
[above](#2-computing-the-total-binding-affinity))

4. (Optional) List of samples ID

**Command**

```
Rscript compute_delta_tba.R \
    --bed regreg/chr7/regreg.bed \                       # The regulatory regions
    --vcf vcf/expr_dataset/chr7/vcf.gz \                 # The VCF file
    --tba tba/expr_dataset/chr7/tba.tsv \                # The total binding affinity
    --snps tba/expr_dataset/chr7/vcfrider_snps.tsv \     # SNPs file (output of vcf\_rider)
    --threads 1 \                                        # Number of cores for parallelisation
    --output delta_tba/expr_dataset/chr7/delta_tba.tsv \ # Output file
```

### 2. Performing the summary-level TWAS

The goal is to find a chain of correlations:
```
SNP->TBA->Expression->Phenotype
```

**Input**

4. BED file with the regulatory regions (obtained previously)

1. The delta-TBA (obtained in the previous step)

3. The weights of the affinities on gene expression (obtained in the
training)

2. Summary statistics (a GWAS Z-score for each SNP)

5. A plink-BED file from a reference population (used to compute LD)

**Remarks**

For best results, we recommend to impute and munge the summary
statistics with [fizi](https://github.com/bogdanlab/fizi). The summary
statistics themselves can be obtained in several ways: they can be
downloaded from the web or computed in house by performing a GWAS with
plink (or any other tool to perform GWA studies).

The summary-level TWAS also requires information about the linkage
disequilibrium (LD) between SNPs. This information is passed to TReX
through a plink-BED file from a reference population (note that there
is a difference between UCSC-BED files, which store coordinates, and
plink-BED files, which store genotypes). For example, you can download
the plink-BED of the 1000 Genomes project, or convert the VCF of the
training data set into plink-BED format.

**Command**

```
Rscript summary_twas.R \
    --regreg regreg/chr7/trex_regreg.bed \              # Regulatory regions
    --tba delta_tba/expr_dataset/chr7/delta_tba.Rds \   # Delta-TBA
    --weights training/chr7/fit/trex_weights.tsv \      # Weights
    --zscores sumstats/phenotype/sumstats.tsv \         # Summary statistics
    --ld genotypes/ld_dataset/chr7.bed \                # Reference LD
    --threads 1 \                                       # Number of cores for parallilisation
    --output summary_twas/phenotype/trex_stwas.tsv      # Output file
```

# References

[manuscript in preparation]

# Further help

Please feel free to open an issue in the [GitHub
repository](https://github.com/fmarotta/TReX) or contact the authors by
email for further help.
