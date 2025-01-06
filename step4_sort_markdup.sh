#!/usr/bin/env bash
###############################################################################
# step4_sort_markdup.sh
# STEP 4: Sort & Mark Duplicates
###############################################################################

echo "=== STEP 4: Sort & Mark Duplicates ==="

mkdir -p "${OUTPUT_DIR}/${BAM_DIR}"

for bam_file in "${OUTPUT_DIR}/${ALIGN_DIR}"/*.bam
do
  sample_name=$(basename "$bam_file" .bam)

  # Sort
  "${SAMTOOLS}" sort -@ "${THREADS}" -o "${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_sorted.bam" "${bam_file}"

  # Mark duplicates (using Picard)
  java -jar "${PICARD}" MarkDuplicates \
    I="${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_sorted.bam" \
    O="${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_sorted_markdup.bam" \
    M="${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_markdup_metrics.txt" \
    CREATE_INDEX=true

  # Remove intermediate file
  rm "${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_sorted.bam"
done
