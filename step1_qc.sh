#!/usr/bin/env bash
###############################################################################
# step1_qc.sh
# STEP 1: Quality Control with FastQC and MultiQC
###############################################################################

echo "=== STEP 1: Quality Control ==="

# Create QC directory if not exists
mkdir -p "${OUTPUT_DIR}/${QC_DIR}"

for fq in "${FASTQ_DIR}"/*.fastq.gz
do
  "${FASTQC}" "${fq}" --threads "${THREADS}" --outdir "${OUTPUT_DIR}/${QC_DIR}"
done

# Optionally run MultiQC for consolidated reports
"${MULTIQC}" "${OUTPUT_DIR}/${QC_DIR}" -o "${OUTPUT_DIR}/${QC_DIR}"
