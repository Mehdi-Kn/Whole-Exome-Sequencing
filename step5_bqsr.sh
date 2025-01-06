#!/usr/bin/env bash
###############################################################################
# step5_bqsr.sh
# STEP 5: Base Quality Score Recalibration (BQSR)
###############################################################################

echo "=== STEP 5: Base Quality Score Recalibration (BQSR) ==="

for bam_file in "${OUTPUT_DIR}/${BAM_DIR}"/*_sorted_markdup.bam
do
  sample_name=$(basename "$bam_file" .bam)

  # Create recalibration table
  "${GATK}" BaseRecalibrator \
    -R "${REFERENCE}" \
    -I "${bam_file}" \
    --known-sites "${DBSNP}" \
    --known-sites "${KNOWN_INDELS}" \
    -O "${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_recal_data.table"

  # Apply recalibration
  "${GATK}" ApplyBQSR \
    -R "${REFERENCE}" \
    -I "${bam_file}" \
    --bqsr-recal-file "${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_recal_data.table" \
    -O "${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_final.bam"

  # Index final BAM
  "${SAMTOOLS}" index "${OUTPUT_DIR}/${BAM_DIR}/${sample_name}_final.bam"
done
