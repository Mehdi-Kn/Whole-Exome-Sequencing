#!/usr/bin/env bash
###############################################################################
# step2_trimming.sh
# STEP 2: Adapter Trimming with Trim Galore
###############################################################################

echo "=== STEP 2: Adapter Trimming ==="

mkdir -p "${OUTPUT_DIR}/${TRIM_DIR}"

for fq1 in "${FASTQ_DIR}"/*_R1*.fastq.gz
do
  # Find matching R2 file based on your naming convention
  fq2="${fq1/_R1/_R2}"
  sample_name=$(basename "$fq1" | sed 's/_R1.*//')

  "${TRIMGALORE}" \
    --paired \
    --quality 20 \
    --cores "${THREADS}" \
    --output_dir "${OUTPUT_DIR}/${TRIM_DIR}" \
    "$fq1" "$fq2"
done
