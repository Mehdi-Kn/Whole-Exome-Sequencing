#!/usr/bin/env bash
###############################################################################
# step6_variant_calling.sh
# STEP 6: Variant Calling (GATK HaplotypeCaller)
###############################################################################

echo "=== STEP 6: Variant Calling with HaplotypeCaller ==="

mkdir -p "${OUTPUT_DIR}/${VCF_DIR}"

for final_bam in "${OUTPUT_DIR}/${BAM_DIR}"/*_final.bam
do
  sample_name=$(basename "$final_bam" _final.bam)

  "${GATK}" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${final_bam}" \
    -O "${OUTPUT_DIR}/${VCF_DIR}/${sample_name}.g.vcf.gz" \
    -ERC GVCF \
    -L "${TARGET_BED}"   # optional but recommended if you have an exome bed
done
