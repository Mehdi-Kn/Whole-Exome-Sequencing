#!/usr/bin/env bash
###############################################################################
# step8_annotation.sh
# STEP 8: Annotation (e.g., ANNOVAR)
###############################################################################

echo "=== STEP 8: Annotation ==="

mkdir -p "${OUTPUT_DIR}/${ANN_DIR}"

for vcf_file in "${OUTPUT_DIR}/${VCF_DIR}"/*_filtered.vcf.gz
do
  sample_name=$(basename "$vcf_file" _filtered.vcf.gz)

  # Decompress into .vcf for ANNOVAR
  "${BCFTOOLS}" view "${vcf_file}" > "${OUTPUT_DIR}/${ANN_DIR}/${sample_name}_filtered.vcf"

  # Convert VCF to ANNOVAR input
  perl "${ANNOVAR_DIR}/convert2annovar.pl" \
    -format vcf4 \
    "${OUTPUT_DIR}/${ANN_DIR}/${sample_name}_filtered.vcf" \
    > "${OUTPUT_DIR}/${ANN_DIR}/${sample_name}.avinput"

  # Run table_annovar.pl with your chosen databases
  perl "${ANNOVAR_DIR}/table_annovar.pl" \
    "${OUTPUT_DIR}/${ANN_DIR}/${sample_name}.avinput" \
    "${ANNOVAR_DIR}/humandb/" \
    -buildver hg38 \
    -out "${OUTPUT_DIR}/${ANN_DIR}/${sample_name}" \
    -remove \
    -protocol refGene,cytoBand,gnomad211_exome \
    -operation g,r,f \
    -nastring . \
    -vcfinput
done
