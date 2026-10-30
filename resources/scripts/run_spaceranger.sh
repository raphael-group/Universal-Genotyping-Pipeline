#!/bin/bash
# Run Space Ranger count on a Visium HD 3' slide with BAM output, the input UGP genotypes from
#
# Runpeng Luo (2026-10-29)
#
# Dependencies:
#   Space Ranger >= 4.0 (Visium HD 3' needs 4.0+), a Space Ranger transcriptome (e.g. refdata-gex-GRCh38-2024-A)
#
# Usage:
#   bash resources/scripts/run_spaceranger.sh <spaceranger> <transcriptome> <id> <slide-area> \
#           <fastqs> <cytaimage> <image> <slidefile> <output_dir>
#           spaceranger    # spaceranger binary
#           transcriptome  # Space Ranger reference dir; its fasta/genes match UGP reference/gtf_file
#           id             # Space Ranger run id
#           slide-area     # slide serial and capture area, e.g. H1-CWRFYX3-A1
#           fastqs         # dir of the unpacked fastqs
#           cytaimage      # CytAssist image (.tif)
#           image          # H&E tissue image (.btf)
#           slidefile      # slide layout file (.vlf)
#           output_dir     # Space Ranger output dir; results under <output_dir>/outs
#
# Inputs
#   fastqs, cytaimage, image, slidefile: the slide's Space Ranger inputs
# Parameters
#   LOCALCORES, LOCALMEM: environment overrides of cores and GB of memory (default 48, 64)
# Outputs:
#   <output_dir>/outs/possorted_genome_bam.bam{,.bai}: 2um-barcode BAM (CB = square_002um)
#   <output_dir>/outs/barcode_mappings.parquet: square_002um -> square_008um/square_016um/cell
#   <output_dir>/outs/binned_outputs/square_{002,008,016}um/: filtered_feature_bc_matrix{,.h5}, spatial/
#   <output_dir>/outs/spatial/: slide images
# Notes/References:
#   --create-bam=true is required: UGP piles up allele counts from the BAM.
#   Space Ranger count: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/space-ranger-count
#   Outputs: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview
set -euo pipefail

if [ $# -ne 9 ]; then
    sed -n 2,21p "$0"
    exit 1
fi
spaceranger=$1
transcriptome=$2
id=$3
slide_area=$4
fastqs=$5
cytaimage=$6
image=$7
slidefile=$8
output_dir=$9

LOCALCORES=${LOCALCORES:-48}
LOCALMEM=${LOCALMEM:-64}

echo "[$(date '+%F %T')] spaceranger count ${id} (${slide_area}) -> ${output_dir}"
"$spaceranger" count \
    --id="$id" \
    --slide="${slide_area%-*}" \
    --area="${slide_area##*-}" \
    --transcriptome="$transcriptome" \
    --fastqs="$fastqs" \
    --cytaimage="$cytaimage" \
    --image="$image" \
    --slidefile="$slidefile" \
    --create-bam=true \
    --output-dir="$output_dir" \
    --localcores="$LOCALCORES" \
    --localmem="$LOCALMEM"
