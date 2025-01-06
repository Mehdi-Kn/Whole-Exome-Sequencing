#!/usr/bin/env bash
###############################################################################
# step3_alignment.sh
# STEP 3: Alignment with BWA MEM
###############################################################################

echo "=== STEP 3: Alignment with BWA ==="

mkdir -p "${OUTPUT_DIR}/${ALIGN_DIR}"

for fq1_trimmed in "${OUTPUT_DIR}/${TRIM_DIR}"/*_R1_*.fq.gz
do
  fq2_trimmed="${fq1_trimmed/_R1_/_R2_}"
  sample_name=$(basename "$fq1_trimmed" | sed 's/_R1_.*//')

  "${BWA}" mem -t "${THREADS}" \
    "${REFERENCE}" \
    "${fq1_trimmed}" \
    "${fq2_trimmed}" \
  | "${SAMTOOLS}" view -bS - \
  > "${OUTPUT_DIR}/${ALIGN_DIR}/${sample_name}.bam"
done
