# Whole-Exome-Sequencing Pipeline

Welcome to my **Whole Exome Sequencing (WES) Pipeline**! 
- This pipeline is designed to take you from raw exome sequencing data all the way to meaningful variant discovery.
- It is streamlined, easy to use, and follows best practices (e.g., GATK guidelines) for high-quality, reproducible results.

## Table of Contents
1. [Biological Overview](#biological-overview)  
2. [Pipeline Steps](#pipeline-steps)  
   1. [Quality Control (QC)](#1-quality-control-qc)  
   2. [Adapter Trimming](#2-adapter-trimming)  
   3. [Alignment/Mapping](#3-alignmentmapping)  
   4. [Marking Duplicates](#4-marking-duplicates)  
   5. [Base Quality Score Recalibration (BQSR)](#5-base-quality-score-recalibration-bqsr)  
   6. [Variant Calling](#6-variant-calling)  
   7. [Variant Filtering](#7-variant-filtering)  
   8. [Annotation](#8-annotation)  
3. [Usage Instructions](#usage-instructions)  
4. [Dependencies](#dependencies)  
5. [License](#license)  
6. [Contact](#contact)

---

## Biological Overview

Whole Exome Sequencing (WES) targets the **protein-coding regions** of the genome (exons), where the majority of known disease-causing variants reside. By focusing on these regions, WES provides a cost-effective approach to detect both common and rare variants linked to human diseases or traits.

In this pipeline, you’ll move from **raw FASTQ files** to **high-confidence variant calls (VCFs)**, ready for downstream analyses or clinical interpretation.

---

## Pipeline Steps

### 1. Quality Control (QC)
**Biological Meaning**  
Ensuring that your raw reads are of high quality is crucial for accurate downstream analyses. Issues like low-quality bases or adapter contamination can obscure true variants.

**Overview**  
Common QC tools (e.g., **FastQC**, **MultiQC**) assess base quality scores, GC content, adapter contamination, and other metrics. Reviewing QC reports early helps detect problems before they propagate through the pipeline.

**Pro Tip**  
If you notice excessive adapter content or poor average quality in the last few cycles, you may need more aggressive trimming or a new library prep strategy.

---

### 2. Adapter Trimming
**Biological Meaning**  
Removing adapter sequences and low-quality bases helps improve alignment rates and overall variant calling accuracy.

**Overview**  
Tools like **Trim Galore** or **Trimmomatic** detect and remove adapter sequences and low-quality bases. This ensures subsequent mapping focuses on the highest-quality, most informative parts of each read.

**Pro Tip**  
Balance is key. Over-trimming can result in very short reads (losing coverage), while under-trimming leaves noisy bases that can lead to false variant calls.

---

### 3. Alignment/Mapping
**Biological Meaning**  
Accurate mapping of reads to a reference genome (e.g., **hg19** or **hg38**) is foundational for reliable variant detection. Misalignments can introduce false positives or mask true variants.

**Overview**  
Common tools (e.g., **BWA**, **Bowtie2**) align reads to the reference. This process compares each read against the reference sequence, creating a **BAM** file that stores read positions and quality data.

**Pro Tip**  
Use the same reference genome consistently (e.g., if you’re using **hg38**, all annotation files, BED files, and known variants should also be **hg38**).

---

### 4. Marking Duplicates
**Biological Meaning**  
During library preparation, **PCR amplification** can create duplicate reads that don’t represent unique DNA fragments. Marking duplicates prevents overcounting in variant calling.

**Overview**  
Tools like **Picard MarkDuplicates** or **GATK MarkDuplicatesSpark** identify reads with identical alignment coordinates. These are flagged in the BAM file so they’re not double-counted.

**Pro Tip**  
A very high duplication rate may indicate low library complexity. You might need to optimize your library prep method or input DNA amount.

---

### 5. Base Quality Score Recalibration (BQSR)
**Biological Meaning**  
Sequencing machines can introduce systematic errors in base quality scores, leading to false variant calls. **BQSR** aligns quality scores with empirically observed error rates.

**Overview**  
Using **GATK BaseRecalibrator** and **ApplyBQSR**, you’ll correct base quality scores based on known sites (e.g., **dbSNP**). This improves overall variant call accuracy.

**Pro Tip**  
Ensure the known variant databases (dbSNP, Mills, 1000G, etc.) match your reference genome build (hg19 vs. hg38). Mismatched references can lead to incorrect recalibrations.

---

### 6. Variant Calling
**Biological Meaning**  
Variant calling identifies genomic alterations (SNPs, indels, etc.) within the exons. These variants may be implicated in disease or other phenotypic traits.

**Overview**  
Tools like **GATK HaplotypeCaller** or **FreeBayes** examine read evidence at each position to determine if a variant exists. For WES, calling can be restricted to exome capture regions for speed and efficiency.

**Pro Tip**  
For multi-sample projects, consider using **joint calling** approaches (e.g., **GenotypeGVCFs** in GATK) to improve accuracy across a cohort.

---

### 7. Variant Filtering
**Biological Meaning**  
Filtering out low-quality variant calls refines the candidate list for downstream analysis. This step helps distinguish **true positives** from artifacts.

**Overview**  
You can set thresholds (e.g., **QD**, **FS**, **MQ**) using **GATK VariantFiltration** or other tools. Well-chosen filters significantly reduce false positives.

**Pro Tip**  
Don’t be too aggressive—over-filtering can remove real variants. Base your thresholds on read depth, coverage, and prior knowledge of your samples.

---

### 8. Annotation
**Biological Meaning**  
Annotation interprets each variant’s functional impact (e.g., synonymous, nonsynonymous, frameshift) and checks databases like **gnomAD**, **ClinVar**, or **COSMIC**.

**Overview**  
Tools like **ANNOVAR**, **VEP**, or **SnpEff** map variants to genes and annotate them with data from public databases. This provides insight into the clinical or biological significance of variants.

**Pro Tip**  
Regularly update your annotation databases. New findings may alter the classification of existing variants.

---

## Usage Instructions

1. **Clone the Repository**  
   ```bash
   git clone https://github.com/Mehdi-Kn/Whole-Exome-Sequencing.git
   cd Whole-Exome-Sequencing
   ```

2. **Prepare Your Data**  
   - Place raw FASTQ files in a `data/raw/` folder (or similar).  
   - Ensure you have a reference genome (e.g., `hg38.fa`) with a corresponding `.fai`, `.dict`, and BWA index in a `reference/` folder.  
   - Download known variant databases (dbSNP, Mills, 1000G, etc.) for BQSR if needed.

3. **Install or Load Tools**  
   - **FastQC**, **MultiQC**  
   - **Trim Galore** (or **Trimmomatic**)  
   - **BWA**, **SAMtools**  
   - **Picard Tools**  
   - **GATK** (>= 4.x)  
   - **bcftools** (for VCF manipulation)  
   - **ANNOVAR** or **VEP** (for annotation)

4. **Run Each Step** (Example Commands)  
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
   #   e.g., picard MarkDuplicates or GATK MarkDuplicatesSpark
   
   # 5. BQSR
   #   e.g., gatk BaseRecalibrator / ApplyBQSR
   
   # 6. Variant Calling
   gatk HaplotypeCaller -R reference/hg38.fa -I sample_final.bam -O sample.g.vcf.gz ...
   
   # 7. Variant Filtering
   #   e.g., gatk VariantFiltration ...
   
   # 8. Annotation
   #   e.g., ANNOVAR or VEP...
   ```

   Adjust paths, filenames, and parameters as needed for your environment.

5. **Review Results**  
   - **QC Reports**: `results/qc/`  
   - **Trimmed FASTQs**: `data/trimmed/`  
   - **Aligned BAMs**: `results/aligned/`  
   - **Raw VCFs**: `results/variants/`  
   - **Filtered or Annotated VCFs**: `results/filtered/`, `results/annotated/`  

---

## Dependencies

- **Operating System**: Linux (Ubuntu, CentOS, etc.)  
- **Bioinformatics Tools**:  
  - **fastqc**, **multiqc**, **trim_galore** (or **trimmomatic**), **bwa**, **samtools**  
  - **picard**, **gatk** (>= 4.x)  
  - **bcftools**  
  - **annovar** or **vep** (for annotation)  
- **Reference Files**:  
  - Genome FASTA (hg19 or hg38) + `.fai`, `.dict`, and BWA index  
  - Known sites (dbSNP, Mills, 1000G) for BQSR (matching reference build)  
- **Installation Notes**:  
  - GATK and Picard often require Java 8+.  
  - Ensure reference genome and known sites are consistently the same assembly.

---

## License

This project is licensed under the **MIT License**. You are free to modify or distribute it as needed.

---

## Contact

**Author**: Knidiri Mehdi  
**Email**: m.knidiri70@gmail.com  

Feel free to open an issue or submit a pull request if you have suggestions or improvements!

---

### Final Notes

- **Customization**: Adjust sample names, threshold values, and post-processing steps to match your experimental design.  
- **Capture Kits**: If using an exome capture kit, specify the capture BED (e.g., `-L capture.bed`) in GATK commands to restrict analysis to target regions.  
- **Multi-Sample Workflows**: For cohort studies, consider **joint variant calling** (e.g., GATK `GenomicsDBImport` + `GenotypeGVCFs`).  
- **Downstream Analysis**: After obtaining final VCFs, you can proceed to variant prioritization, functional interpretation, and clinical relevance assessments.

**GOOD LUCK** with your analyses!
