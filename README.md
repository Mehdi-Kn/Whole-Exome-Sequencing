# Whole-Exome-Sequencing Pipeline

Welcome to my **WES Pipeline**! This pipeline is designed to take you from raw exome sequencing data all the way to meaningful variant discovery. 
- It's streamlined, easy to use, and includes tips and tricks to simplify your analysis. 
- Following best practices (e.g., GATK guidelines), this pipeline ensures high-quality, reproducible results.

## Table of Contents

1. [Biological Overview](#biological-overview)  
2. [Pipeline Steps](#pipeline-steps)  
   - [1. Quality Control (QC)](#1-quality-control-qc)  
   - [2. Adapter Trimming](#2-adapter-trimming)  
   - [3. Alignment/Mapping](#3-alignmentmapping)  
   - [4. Marking Duplicates](#4-marking-duplicates)  
   - [5. Base Quality Score Recalibration (BQSR)](#5-base-quality-score-recalibration-bqsr)  
   - [6. Variant Calling](#6-variant-calling)  
   - [7. Variant Filtering](#7-variant-filtering)  
   - [8. Annotation](#8-annotation)  
3. [Usage Instructions](#usage-instructions)  
4. [Dependencies](#dependencies)  
5. [License](#license)  
6. [Contact](#contact)  

---

## Biological Overview

Whole Exome Sequencing (WES) targets the protein-coding regions of the genome (exons), which harbor the majority of known disease-causing variants. By focusing on these regions, WES provides a cost-effective approach to detect both common and rare variants linked to human diseases or traits.

In this pipeline, you’ll move from raw FASTQ files to high-confidence variant calls (VCFs), ready for downstream analyses or clinical interpretation.

---

## Pipeline Steps

### 1. Quality Control (QC)

**Biological Meaning:**  
Ensuring that your raw reads are of high quality is crucial for accurate downstream analyses. Issues like low-quality bases or adapter contamination can obscure true variants.

**Overview:**  
Common QC tools (e.g., **FastQC**, **MultiQC**) assess base quality scores, GC content, adapter contamination, and other metrics. Reviewing QC reports early helps detect problems before they propagate through the pipeline.

**Pro Tip:**  
If you notice excessive adapter content or poor average quality in the last few cycles, you may need more aggressive trimming or a new library prep strategy.

---

### 2. Adapter Trimming

**Biological Meaning:**  
Trimming adapters and removing low-quality bases helps improve alignment rates and overall variant calling accuracy.

**Overview:**  
Tools like **Trim Galore** or **Trimmomatic** can automatically detect and remove adapter sequences and low-quality regions. This ensures that subsequent mapping focuses on the highest-quality, most informative parts of each read.

**Pro Tip:**  
Balance is key. Over-trimming can result in very short reads, losing coverage, while under-trimming leaves noisy bases that can lead to false variant calls.

---

### 3. Alignment/Mapping

**Biological Meaning:**  
Accurate mapping of reads to a reference genome (e.g., hg19 or hg38) is foundational for reliable variant detection. Misalignments can introduce false positives or mask true variants.

**Overview:**  
Common tools (e.g., **BWA**, **Bowtie2**) align reads to the reference. The alignment process compares each read against the reference sequence, creating a **BAM** file that stores read positions and qualities.

**Pro Tip:**  
Make sure your reference genome is consistent throughout the pipeline (e.g., if you’re using hg38, all annotation files, BED files, and known variants should also be hg38).

---

### 4. Marking Duplicates

**Biological Meaning:**  
During library preparation, PCR amplification can create duplicate reads that don’t represent unique DNA fragments. Marking duplicates prevents overcounting.

**Overview:**  
Tools like **Picard MarkDuplicates** or **GATK MarkDuplicatesSpark** identify reads with identical alignment coordinates. These are flagged in the BAM file so they’re not double-counted or misinterpreted during variant calling.

**Pro Tip:**  
High duplication rates can suggest issues with library complexity. If your duplication rate is extremely high, consider adjusting the input DNA amount or library prep method.

---

### 5. Base Quality Score Recalibration (BQSR)

**Biological Meaning:**  
Sequencing machines can introduce systematic errors in base quality scores, which can lead to false variant calls. Recalibrating these scores aligns them with empirical error rates.

**Overview:**  
Using **GATK BaseRecalibrator** and **ApplyBQSR**, you’ll correct base quality scores based on known sites (e.g., dbSNP). Recalibration helps produce more accurate variant calls.

**Pro Tip:**  
Always use the same genome build (e.g., hg19 or hg38) for known sites. Mismatched references can lead to incorrect recalibrations.

---

### 6. Variant Calling

**Biological Meaning:**  
Variant calling identifies genomic alterations—SNPs, indels, etc.—within the exons. This reveals potential mutations driving disease or other phenotypic traits.

**Overview:**  
Tools like **GATK HaplotypeCaller** or **FreeBayes** examine the read evidence at each position to determine if a variant is present. For WES, calling is typically restricted to exome capture regions for improved efficiency.

**Pro Tip:**  
For multi-sample projects, consider performing joint calling or GenotypeGVCFs (if using GATK) for improved accuracy across your cohort.

---

### 7. Variant Filtering

**Biological Meaning:**  
Filtering out low-quality variant calls helps distinguish true positives from artifacts. This step refines the list of candidate variants for downstream analyses.

**Overview:**  
You can set thresholds based on quality metrics (e.g., QD, FS, MQ) using **GATK VariantFiltration** or equivalent tools. Good filtering steps significantly reduce false-positive rates.

**Pro Tip:**  
Don’t be too aggressive with filters—over-filtering can remove true variants. Instead, adjust threshold criteria based on coverage, read depth, and prior knowledge of your samples.

---

### 8. Annotation

**Biological Meaning:**  
Annotation interprets the functional impact of each variant—e.g., whether it’s synonymous, nonsynonymous, frameshift, or located in a regulatory region—and checks databases like gnomAD, ClinVar, or COSMIC.

**Overview:**  
Tools like **ANNOVAR**, **VEP**, or **SnpEff** map your variants to genes and add information from public databases. This helps you understand the clinical or biological significance of variants.

**Pro Tip:**  
Regularly update your annotation databases. Genetic databases evolve quickly, and new information might change the classification of a variant.

---

## Usage Instructions

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/Mehdi-Kn/WES-Pipeline.git
   cd WES-Pipeline
   ```

2. **Prepare Your Data:**
   - Place your raw FASTQ files in a `data/` folder (e.g., `data/raw`).
   - Ensure you have a reference genome (`hg38.fa`) and a corresponding index/dictionary in a `reference/` folder.
   - Download known variant databases (e.g., dbSNP, Mills) if you plan to do BQSR.

3. **Install or Load Tools:**
   - FastQC, MultiQC  
   - Trim Galore or Trimmomatic  
   - BWA, SAMtools  
   - GATK (>= 4.x)  
   - Picard Tools  
   - ANNOVAR or VEP  

4. **Run Each Step:**

   ```bash
   # 1. QC
   fastqc data/raw/*.fastq.gz -o results/qc/
   multiqc results/qc/ -o results/qc/

   # 2. Adapter Trimming
   trim_galore --paired data/raw/sample_R1.fastq.gz data/raw/sample_R2.fastq.gz -o data/trimmed/

   # 3. Alignment
   bwa mem reference/hg38.fa data/trimmed/sample_R1_val_1.fq.gz data/trimmed/sample_R2_val_2.fq.gz \
       | samtools view -bS - > results/aligned/sample.bam

   # 4. Sort & Mark Duplicates
   # ... etc. (Picard or GATK MarkDuplicates)

   # 5. BQSR
   # ... GATK BaseRecalibrator / ApplyBQSR

   # 6. Variant Calling
   gatk HaplotypeCaller -R reference/hg38.fa -I sample_final.bam -O sample.g.vcf.gz --some-params

   # 7. Variant Filtering
   # gatk VariantFiltration ...

   # 8. Annotation
   # ANNOVAR or VEP...
   ```

   Customize these commands to match your file structure and parameters.

5. **Review Results:**
   - **BAM Files:** `results/aligned/`
   - **QC Reports:** `results/qc/`
   - **Raw VCF Files:** `results/variants/`
   - **Filtered/Annotated VCFs:** `results/filtered/`, `results/annotated/`

---

## Dependencies

- **Operating System:** Linux (Ubuntu, CentOS, etc.)  
- **Bioinformatics Tools:**  
  - `fastqc`, `multiqc`, `trim_galore` (or `trimmomatic`), `bwa`, `samtools`  
  - `picard`, `gatk` (>= 4.x)  
  - `bcftools` (for VCF manipulation)  
  - `annovar` or `vep` for annotation  
- **Reference Files:**  
  - `hg19` or `hg38` FASTA + `.fai`, `.dict`, and BWA index  
  - Known sites (dbSNP, Mills, 1000G) for BQSR (matching reference build)

**Installation Notes:**  
- Tools like GATK and Picard often require Java 8+.  
- Ensure your reference genome and known variant VCFs align with the same assembly.  
- For reproducibility, consider using `conda` or Docker containers.

---

## License

This project is licensed under the [MIT License](LICENSE).  
Feel free to modify or distribute as needed.

---

## Contact

For questions, suggestions, or contributions, please reach out to:

**Knidiri Mehdi**  
Email: [m.knidiri70@gmail.com](m.knidiri70@gmail.com)

---

# Final Notes

- **Customization:** Adjust sample names, thresholding, and post-processing steps based on your specific experimental design.  
- **Capture Kits:** If you have a capture kit BED file, pass it to `-L capture.bed` when calling variants to speed up GATK.  
- **Multi-Sample Workflows:** For population or cohort studies, consider joint calling (e.g., GATK GenomicsDBImport and GenotypeGVCFs).  
- **Downstream Analysis:** Once you have final VCFs, explore variant prioritization, functional interpretation, and clinical relevance using specialized databases.  

By following these steps, you’ll be set to run a robust, reproducible WES pipeline. GOOD LUCK .
