#!/bin/bash -x

set -euo pipefail

prefix="${1}"

s3_bam_file="${2}"
bam_filename=$(basename "${s3_bam_file}")
BAM_FOLDER="BAMs"
bam_file="${BAM_FOLDER}/${bam_filename}"

s3_bai_file="${3}"
bai_file="${BAM_FOLDER}/${bam_filename}.bai"

s3_binmap_file="${4}"
binmap_file=$(basename "${s3_binmap_file}")

library_size="${5}"
target_genome="${6}"

if [ ! -f "${binmap_file}" ]; then
    aws s3 cp "${s3_binmap_file}" "${binmap_file}" --quiet
fi

mkdir -p ${BAM_FOLDER}
if [ ! -f "${bam_file}" ]; then
    aws s3 cp "${s3_bam_file}" "${bam_file}" --quiet
    aws s3 cp "${s3_bai_file}" "${bai_file}" --quiet
fi

python3 ninjamap2.py $bam_file $binmap_file $library_size $prefix $target_genome | tee -a ${prefix}.log
