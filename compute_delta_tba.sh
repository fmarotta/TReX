#!/usr/bin/env bash
set -euo pipefail
trap clean INT

# compute_delta_tba.sh
#
# Given reference SNPs, DNA regions, reference sequence and PWM, compute
# the average tba for each allele.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Sep 2020

usage=$(cat <<-END
Usage: $(basename $0) OPTIONS

Compute the Total Binding Affinity (TBA) of each transcription factor for each
regulatory region.

Options:
    -v --vcf=<VCF>              The reference LD dataset as a phased VCF
    -b --bed=<BED_FILE>         Interesting regions
    -r --ref=<FASTA_FILE>       The reference sequence of the genome
    -p --pwm=<TRANSFAC_FILE>    The PWM in transfac format
    -f --bkg=<BKG_FILE>         Background frequencies [default: 'uniform']
    -k --keep=<FILE>            List of samples to keep [default: all of 'em]
    -l --too-large=<FLOAT>      Fraction to discard indels [default: .1]
    -o --outdir=<DIR>           Path to output files [default: ./]
    -j --threads=<N>            Number of cores for parallel computations [default: 1]
    -h --help                   Print this message
END
)

# Default options
fraction=.1
threads=1
outdir="."

# Custom functions
function fvalid()
{
    if [ ! -f "$1" ]; then
        if [ ! -p "$1" ]; then # It is OK if it's a pipe
            echo -e "Error: $1 does not exist.\n" >&2
            echo "$usage" >&2
            exit 1
        fi
    fi
}
function dvalid()
{
    if [ ! -d "$1" ]; then
        echo -e "Info: creating new directory $1.\n" >&2
        if ! mkdir -p $1; then
            echo -e "Error: could not create directory $1.\n"
            exit 1
        fi
    fi
}
function clean()
{
    echo "Removing all output..." >&2
    rm $reduced_vcf $windows $chk $snps $vcfrider_tba $tba
}


# Options string
options=v:b:r:p:f:k:l:j:o:h
longoptions=vcf:,bed:,ref:,pwm:,bkg:,keep:,too-large:,threads:,outdir:,help

# check syntax and acquire options and arguments
PARSER=$(getopt --options=$options --longoptions=$longoptions --name "$0" -- "$@")
eval set -- "$PARSER"

while true; do
    case "$1" in
        -h|--help )
            echo "$usage"
            exit 0
            ;;
        -v|--vcf )
            vcf=$2
            shift 2
            ;;
        -b|--bed )
            bed=$2
            shift 2
            ;;
        -r|--ref )
            ref=$2
            shift 2
            ;;
        -p|--pwm )
            pwm=$2
            shift 2
            ;;
        -f|--bkg )
            bkg=$2
            shift 2
            ;;
        -k|--keep )
            keep=$2
            shift 2
            ;;
        -l|--too-large )
            fraction=$2
            shift 2
            ;;
        -j|--threads )
            threads=$2
            shift 2
            ;;
        -o|--outdir )
            outdir=$2
            shift 2
            ;;
        -- )
            shift
            break
            ;;
        * )
            echo "Unexpected arguments problem." >&2
            exit 1
            ;;
    esac
done


# Check for mandatory arguments
for arg in vcf bed ref pwm outdir; do
    if [ -z ${!arg+x} ]; then
        echo "Error: --$arg is not specified, but it is mandatory." >&2
        exit 1
    # else
    #     echo "$arg is set to ${!arg}"
    fi
done

# Create the uniform background, if necessary
if [ -z ${bkg+x} ]; then
    echo "Info: No background provided, generating uniform one..." >&2
    bkg=$outdir/background.txt
    echo -e "0.25\t0.25\t0.25\t0.25" > $bkg
fi

# Check that files exist
for file in vcf bed ref pwm bkg; do
    fvalid ${!file}
done
dvalid $outdir

# Look for vcf_rider
if which vcf_rider > /dev/null; then
    vcfriderpath=$(dirname $(which vcf_rider))
else
    vcfriderpath=$(find . -name vcf_rider | head -n 1)
    if [ -z "$vcfriderpath" ]; then
        echo "Error: couldn't find vcf_rider. Make sure it is either searchable in" >&2
        echo "the PATH environmental variable or under the current working directory." >&2
        exit 1
    fi
    vcfriderpath=$(dirname $vcfriderpath)
fi

# In-script variables
reduced_vcf="$outdir/filtered_phased_vcf.vcf"
windows="$outdir/delta_windows.bed"
chk="$outdir/delta_indel_stats.tsv"
snps="$outdir/delta_vcfrider_snps.tsv"
vcfrider_tba="$outdir/delta_vcfrider_tba.tsv.gz"
tba="$outdir/delta_tba.Rds"


# Filter the VCF and append a heterozygous individual to make sure that
# each allele is represented
echo "Filtering the VCF..." >&2
if [ -z ${keep+x} ]; then
    vcftools --gzvcf $vcf --bed <(sed 's/chr//' $bed) --keep-INFO-all \
        --recode --stdout \
    | awk -F '\t' -v OFS='\t' '($0 ~ /^##/) {print}
        ($0 ~ /^#CHROM/) {print $0, "PSEUDOIID"}
        ($0 !~ /^#/) {print $0, "0|1"}' > $reduced_vcf
else
    vcftools --gzvcf $vcf --bed <(sed 's/chr//' $bed) --keep-INFO-all \
        --keep $keep --recode --stdout \
    | awk -F '\t' -v OFS='\t' '($0 ~ /^##/) {print}
            ($0 ~ /^#CHROM/) {print $0, "PSEUDOIID"}
            ($0 !~ /^#/) {print $0, "0|1"}' > $reduced_vcf
fi


echo "Finding the windows around the variants..." >&2
# Find max TF motif length
l=$(awk 'BEGIN {MAX = 0}
         ($1 ~ /^[0-9]+$/ && $1 > MAX) {MAX = $1}
         END {print MAX}' $pwm)

# Generate the SNP-based BED
awk -F '\t' -v OFS='\t' -v l=$l \
    'function max(n1, n2) {
        if (n1 > n2)
            return n1
        return n2
    }
    ($0 !~ /^#/) {print $1 ~ /chr.*/ ? $1 : "chr"$1,
                        $2 - l - 1,
                        $2 + (max(length($4), length($5)) - 1) + l,
                        $3}' $reduced_vcf \
| bedtools intersect -a stdin -b $bed \
| uniq > $windows


# Indel stats
echo "Running indel_stats..." >&2
$vcfriderpath/indel_stats $reduced_vcf $windows > $chk


# Compute TBA
echo "Computing the TBA..." >&2
$vcfriderpath/vcf_rider \
  -b <(grep -f <(awk -F '\t' -v OFS='\t' '$6 == "false" {print $1}' $chk) $windows) \
  -p <(/bioinfo/prj/expr_reg_pred/perl/transfac2pcm.pl < $pwm) \
  -f $bkg \
  -v $reduced_vcf \
  -r <(zcat $ref | sed 's/chr//1' | awk -F '\t' -v OFS='\t' '{print toupper($0)}') \
  -a $snps | gzip -c > $vcfrider_tba


# Aggregate the delta TBA
echo "Aggregating the delta TBA..." >&2
Rscript utils/aggregate_delta_tba.R \
    --vcf=$reduced_vcf \
    --tba=$vcfrider_tba \
    --snps=$snps \
    --too-large=$fraction \
    -j $threads \
    -o $tba


# Clear temporary files
# rm $reduced_vcf
# rm $windows
# rm $chk
# rm $snps
# rm $vcfrider_tba
