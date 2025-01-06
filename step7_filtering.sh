#!/usr/bin/env bash
###############################################################################
# step7_filtering.sh
# STEP 7: Variant Filtering
###############################################################################

echo "=== STEP 7: Variant Filtering ==="

for gvcf in "${OUTPUT_DIR}/${VCF_DIR}"/*.g.vcf.gz
do
  sample_name=$(basename "$gvcf" .g.vcf.gz)

  "${GATK}" VariantFiltration \
    -R "${REFERENCE}" \
    -V "${gvcf}" \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    -O "${OUTPUT_DIR}/${VCF_DIR}/${sample_name}_filtered.vcf.gz"
done
