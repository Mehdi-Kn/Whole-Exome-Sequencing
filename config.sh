#!/usr/bin/env bash
###############################################################################
# config.sh
# Define all tool YourOwnPaths, reference files, and default parameters for the WES 
# pipeline. This is sourced by other scripts.
###############################################################################

# -- Tools ---------------------------------------------------------------------
FASTQC="/YourOwnPath/fastqc"
MULTIQC="/YourOwnPath/multiqc"
TRIMGALORE="/YourOwnPath/trim_galore"
BWA="/YourOwnPath/bwa"
SAMTOOLS="/YourOwnPath/samtools"
PICARD="/YourOwnPath/picard.jar"        # or YourOwnPath to MarkDuplicatesSpark if using GATK
GATK="/YourOwnPath/gatk"
BCFTOOLS="/YourOwnPath/bcftools"
ANNOVAR_DIR="/YourOwnPath/annovar"      # or /YourOwnPath/VEP if using VEP

# -- Reference Files & Known Sites --------------------------------------------
REFERENCE="/YourOwnPath/hg38.fa"        # or hg19
DBSNP="/YourOwnPath/dbsnp.vcf.gz"
KNOWN_INDELS="/YourOwnPath/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
TARGET_BED="/YourOwnPath/exome_capture.bed"  # Optional, but recommended for WES

# -- Default Parameters --------------------------------------------------------
THREADS=4  # default threads (can be overridden by command-line args)

# -- Output Directory Structure (relative to --out) ---------------------------
QC_DIR="qc"
TRIM_DIR="trimmed_fastq"
ALIGN_DIR="aligned"
BAM_DIR="bam"
VCF_DIR="vcf"
ANN_DIR="annotation"
