#!/usr/bin/env bash
###############################################################################
# wes_pipeline.sh
# A streamlined Whole Exome Sequencing (WES) pipeline from FASTQ to annotated VCF.
# Written in bash with inline "Tips & Tricks".
#
# Author: Mehdi
# Usage Example:
#   bash wes_pipeline.sh --fastq /path/to/fastq --out /path/to/output --threads 8
#
###############################################################################

# ------------------------------------------------------------------------------
# 1. Parse Command-Line Arguments
# ------------------------------------------------------------------------------
usage() {
  echo "Usage: $0 --fastq <FASTQ_DIR> --out <OUTPUT_DIR> [--threads N]"
  exit 1
}

THREADS=4  # Default number of threads
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --fastq)
      FASTQ_DIR="$2"
      shift; shift
      ;;
    --out)
      OUTPUT_DIR="$2"
      shift; shift
      ;;
    --threads)
      THREADS="$2"
      shift; shift
      ;;
    *)
      usage
      ;;
  esac
done

if [[ -z "${FASTQ_DIR}" || -z "${OUTPUT_DIR}" ]]; then
  usage
fi

###############################################################################
# TIP: Make sure the following paths point to the correct executables or modules
###############################################################################

# ------------------------------------------------------------------------------
# 2. Define Paths to Tools and Reference Files
# ------------------------------------------------------------------------------
# Tools
FASTQC="/path/to/fastqc"
MULTIQC="/path/to/multiqc"
TRIMGALORE="/path/to/trim_galore"
BWA="/path/to/bwa"
SAMTOOLS="/path/to/samtools"
PICARD="/path/to/picard.jar"         # or MarkDuplicatesSpark from GATK
GATK="/path/to/gatk"                 # GATK 4.x
BCFTOOLS="/path/to/bcftools"
ANNOVAR_DIR="/path/to/annovar"       # or /path/to/VEP if you use VEP

# Reference & Known Sites
REFERENCE="/path/to/hg38.fa"         # e.g., hg38 or hg19
DBSNP="/path/to/dbsnp.vcf.gz"        # dbSNP (matching reference build)
KNOWN_INDELS="/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
TARGET_BED="/path/to/exome_capture.bed"  # Exome capture kit BED (optional but recommended)

# ------------------------------------------------------------------------------
# 3. Create Output Directories
# ------------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}/qc"
mkdir -p "${OUTPUT_DIR}/trimmed_fastq"
mkdir -p "${OUTPUT_DIR}/aligned"
mkdir -p "${OUTPUT_DIR}/bam"
mkdir -p "${OUTPUT_DIR}/vcf"
mkdir -p "${OUTPUT_DIR}/annotation"

###############################################################################
# TIP: Make sure your reference is indexed:
#   - bwa index /path/to/hg38.fa
#   - samtools faidx /path/to/hg38.fa
#   - gatk CreateSequenceDictionary -R /path/to/hg38.fa -O /path/to/hg38.dict
# And that all known site VCFs are bgzipped and indexed by tabix.
###############################################################################

# ------------------------------------------------------------------------------
# 4. STEP 1: Quality Control (QC)
# ------------------------------------------------------------------------------
echo "=== STEP 1: Quality Control ==="
for fq in "${FASTQ_DIR}"/*.fastq.gz
do
  "${FASTQC}" "$fq" --threads "${THREADS}" --outdir "${OUTPUT_DIR}/qc"
done

# Optionally run MultiQC for consolidated reports
"${MULTIQC}" "${OUTPUT_DIR}/qc" -o "${OUTPUT_DIR}/qc"

###############################################################################
# TIP: Examine the MultiQC HTML report for adapter content, poor quality tails, or
# unusual GC distributions. Adjust trimming parameters as needed.
###############################################################################

# ------------------------------------------------------------------------------
# 5. STEP 2: Adapter Trimming
# ------------------------------------------------------------------------------
echo "=== STEP 2: Adapter Trimming ==="
for fq1 in "${FASTQ_DIR}"/*_R1*.fastq.gz
do
  # Identify the matching R2 file (assuming naming convention: *_R1* and *_R2*)
  fq2="${fq1/_R1/_R2}"
  sample_name=$(basename "$fq1" | sed 's/_R1.*//')
  
  # Using Trim Galore (supports Cutadapt under the hood)
  "${TRIMGALORE}" \
    --paired \
    --quality 20 \
    --cores "${THREADS}" \
    --output_dir "${OUTPUT_DIR}/trimmed_fastq" \
    "$fq1" "$fq2"
done

###############################################################################
# TIP: If you have single-end reads, remove the --paired option.
# Adjust the --quality threshold or other trimming params as needed.
###############################################################################

# ------------------------------------------------------------------------------
# 6. STEP 3: Alignment (BWA MEM)
# ------------------------------------------------------------------------------
echo "=== STEP 3: Alignment with BWA ==="
for fq1_trimmed in "${OUTPUT_DIR}/trimmed_fastq"/*_R1_*.fq.gz
do
  fq2_trimmed="${fq1_trimmed/_R1_/_R2_}"
  sample_name=$(basename "$fq1_trimmed" | sed 's/_R1_.*//')
  
  # Align reads, convert to BAM on the fly
  "${BWA}" mem -t "${THREADS}" \
    "${REFERENCE}" \
    "$fq1_trimmed" \
    "$fq2_trimmed" \
  | "${SAMTOOLS}" view -bS - \
  > "${OUTPUT_DIR}/aligned/${sample_name}.bam"
done

###############################################################################
# TIP: Check alignment statistics (e.g., samtools flagstat) to confirm high
# alignment rates. If rates are low, ensure your reference or read orientation
# is correct.
###############################################################################

# ------------------------------------------------------------------------------
# 7. STEP 4: Sort & Mark Duplicates
# ------------------------------------------------------------------------------
echo "=== STEP 4: Sort & Mark Duplicates ==="
for bam_file in "${OUTPUT_DIR}/aligned"/*.bam
do
  sample_name=$(basename "$bam_file" .bam)
  
  # Sort
  "${SAMTOOLS}" sort -@ "${THREADS}" -o "${OUTPUT_DIR}/bam/${sample_name}_sorted.bam" "$bam_file"
  
  # Mark duplicates (using Picard)
  java -jar "${PICARD}" MarkDuplicates \
    I="${OUTPUT_DIR}/bam/${sample_name}_sorted.bam" \
    O="${OUTPUT_DIR}/bam/${sample_name}_sorted_markdup.bam" \
    M="${OUTPUT_DIR}/bam/${sample_name}_markdup_metrics.txt" \
    CREATE_INDEX=true

  # Remove intermediate file
  rm "${OUTPUT_DIR}/bam/${sample_name}_sorted.bam"
done

###############################################################################
# TIP: If duplicates are extremely high (>50%), you may have library complexity
# or PCR issues. Check duplication metrics in the Picard output.
###############################################################################

# ------------------------------------------------------------------------------
# 8. STEP 5: Base Quality Score Recalibration (BQSR)
# ------------------------------------------------------------------------------
echo "=== STEP 5: BaseRecalibrator & ApplyBQSR ==="
for bam_file in "${OUTPUT_DIR}/bam"/*_sorted_markdup.bam
do
  sample_name=$(basename "$bam_file" .bam)

  # Create recalibration table
  "${GATK}" BaseRecalibrator \
    -R "${REFERENCE}" \
    -I "${bam_file}" \
    --known-sites "${DBSNP}" \
    --known-sites "${KNOWN_INDELS}" \
    -O "${OUTPUT_DIR}/bam/${sample_name}_recal_data.table"

  # Apply recalibration
  "${GATK}" ApplyBQSR \
    -R "${REFERENCE}" \
    -I "${bam_file}" \
    --bqsr-recal-file "${OUTPUT_DIR}/bam/${sample_name}_recal_data.table" \
    -O "${OUTPUT_DIR}/bam/${sample_name}_final.bam"

  # Index final BAM
  "${SAMTOOLS}" index "${OUTPUT_DIR}/bam/${sample_name}_final.bam"
done

###############################################################################
# TIP: BQSR corrects systematic biases in base quality scores. The known-sites
# VCF must match your reference build (hg19 vs. hg38).
###############################################################################

# ------------------------------------------------------------------------------
# 9. STEP 6: Variant Calling (GATK HaplotypeCaller)
# ------------------------------------------------------------------------------
echo "=== STEP 6: Variant Calling with HaplotypeCaller ==="
for final_bam in "${OUTPUT_DIR}/bam"/*_final.bam
do
  sample_name=$(basename "$final_bam" _final.bam)

  "${GATK}" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${final_bam}" \
    -O "${OUTPUT_DIR}/vcf/${sample_name}.g.vcf.gz" \
    -ERC GVCF \
    -L "${TARGET_BED}"   # Optional, recommended if you have an exome capture bed
done

###############################################################################
# TIP: For multiple samples, use GVCF mode + 'GenotypeGVCFs' to combine them:
#   gatk GenotypeGVCFs -R hg38.fa -V sample1.g.vcf.gz -V sample2.g.vcf.gz ... \
#       -O cohort_combined.vcf.gz
###############################################################################

# ------------------------------------------------------------------------------
# 10. STEP 7: Variant Filtering
# ------------------------------------------------------------------------------
echo "=== STEP 7: Variant Filtering ==="
for gvcf in "${OUTPUT_DIR}/vcf"/*.g.vcf.gz
do
  sample_name=$(basename "$gvcf" .g.vcf.gz)

  "${GATK}" VariantFiltration \
    -R "${REFERENCE}" \
    -V "${gvcf}" \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    -O "${OUTPUT_DIR}/vcf/${sample_name}_filtered.vcf.gz"
done

###############################################################################
# TIP: Adjust filter expressions to match GATK best practices or your project's
# requirements (e.g., QD < 2.0, FS > 60, MQ < 40).
###############################################################################

# ------------------------------------------------------------------------------
# 11. STEP 8: Annotation (e.g., ANNOVAR)
# ----------------------------
