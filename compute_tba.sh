#!/usr/bin/env bash
set -euo pipefail
trap clean INT

# compute_tba.sh
#
# Given bed, vcf and pwm, compute the tba.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Feb 2020

usage=$(cat <<-END
Usage: $(basename $0) OPTIONS

Compute the Total Binding Affinity (TBA) of each transcription factor for each
regulatory region.

Options:
    -v --vcf=<VCF_FILE>
    -b --bed=<BED_FILE>
    -r --ref=<FASTA_FILE>
    -p --pwm=<TRANSFAC_FILE>
    -f --bkg=<BKG_FILE>         Background frequencies [default: 'uniform']
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
    rm $chk $snps $vcfrider_tba $sorted_tba $tba
}


# Options string
options=v:b:r:p:f:l:j:o:h
longoptions=vcf:,bed:,ref:,pwm:,bkg:,too-large:,threads:,outdir:,help

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
if [ -z ${bkg+x} ]; then echo "Info: No background provided, generating uniform one..." >&2
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

# Look for blort
if which blort > /dev/null; then
    blortpath=$(dirname $(which blort))
else
    blortpath=$(find . -name blort | head -n 1)
    if [ -z "$blortpath" ]; then
        echo "Error: couldn't find blort. Make sure it is either searchable in" >&2
        echo "the PATH environmental variable or under the current working directory." >&2
        exit 1
    fi
    blortpath=$(dirname $blortpath)
fi

# In-script variables
chk="$outdir/indel_stats.tsv"
snps="$outdir/vcfrider_snps.tsv"
vcfrider_tba="$outdir/vcfrider_tba.tsv"
sorted_tba="$outdir/sorted_tba.tsv"
tba="$outdir/tba.tsv"


# Indel stats
echo "Running indel_stats..." >&2
$vcfriderpath/indel_stats $vcf $bed > $chk


# Compute TBA
echo "Computing the TBA..." >&2
$vcfriderpath/vcf_rider \
  -b <(grep -f <(awk -F '\t' -v OFS='\t' '$6 == "false" {print $1}' $chk) $bed) \
  -p <(perl utils/transfac2pcm.pl < $pwm) \
  -f $bkg \
  -v $vcf \
  -r <(zcat $ref | sed 's/chr//1' | awk -F '\t' -v OFS='\t' '{print toupper($0)}') \
  -a $snps > $vcfrider_tba


# Associate regions
echo "Associating the regions..." >&2
awk -F '\t' -v OFS='\t' '{split($4, a, "_"); print a[1], $4}' $bed | \
sort -t $'\t' -k2,2 | \
join -t $'\t' -1 2 -2 1 -o 1.1,1.2,2.2,2.3,2.4,2.5,2.6,2.7 - \
  <($blortpath/blort $vcfrider_tba) \
> $sorted_tba


# Aggregate TBA
echo "Aggregating the TBA..." >&2
Rscript utils/aggregate_tba.R --too-large=$fraction \
  -j $threads $sorted_tba $snps > $tba


# Clear temporary files
rm $vcfrider_tba
rm $sorted_tba

