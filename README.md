# Bulk RNA-Seq Analysis

🎯 Project Overview
This project demonstrates a complete Bulk RNA-Seq workflow, from raw sequencing data up to differential expression analysis using DESeq2. We analyzed the effect of  **LCOR overexpression** and how it affects gene transcription in two cell lines that differ in nuclear receptor status: 
- **MCF7 (NR-positive)** — nuclear receptor–positive
- **MDA-MB-231 (NR-negative)** — nuclear receptor–negative  
LCOR acts both as a transcriptional corepressor and activator, and this dataset enables comparison of its regulatory activity in nuclear receptor–dependent versus independent contexts.

For each cell line we checked overexpression and normal behaviour of LCOR protein

📂 Dataset Information:
- **GEO Accession:** GSE292767
- **Experiment Type:** Expression profiling by high-throughput sequencing
- **Description:** RNA-Seq experiment in breast cancer cell lines to check the effect of LCOR overexpression in MCF7 and MDA-MB-231.



## Quick Navigation

- [1. Downloading necessary tools](#1-downloading-necessary-tools)
- [2. Fetching SRA Files](#2-fetching-sra-files)
- [3. FastQC - Quality Control](#3-fastqc-quality-control)
- [4. Trimming Reads](#4-trimming-reads)
- [5. Post-trimming Quality Control](#5-post-trimming-quality-control)
- [6. Reference Genome Preparation (Indexing)](#6-reference-genome-preparation-indexing)
- [7. Alignment/Mapping](#7-alignment-mapping)
- [8. SAM/BAM Processing (Sorting/Indexing)](#8-sam-bam-processing-sorting-indexing)
- [9. Assessing Alignment Quality](#9-assessing-alignment-quality)
- [10. Feature Counting (Read Quantification)](#10-feature-counting-read-quantification)
- [11. Downstream Analysis](#11-downstream-analysis)

---

## Workflow

---

## 1. Downloading necessary tools

---

## 2. Fetching SRA Files
- Download SRA files using the SRA toolkit.
- Convert the SRA files into FASTQ files using fastq-dump.
- Save the file in .gz format i.e. gzip file so that it would take less space.

```bash
#Download SRA files
prefetch SRR32858437 SRR32858438 SRR32858439 SRR32858440

#Covert SRA files to FASTQ files
fastq-dump --outdir FASTQ_files --gzip --skip-technical-reads \
--readids --read-filter pass --dumpbase --split-3 --clip \
SRR32858437/SRR32858437.sra
```
**Tools:** `sra-tools` (`prefetch`, `fastq-dump`)

---

## 3. FastQC - Quality Control
- Check the quality of raw sequencing reads using FastQC.

```bash
#Make a directory to save results
mkdir -p QC_results

#Run FASTQC
fastqc FASTQ_files/*.fastq.gz -o QC_results/ --threads 8
```
**Tool:** [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

---

## 4. Trimming Reads (optional)
- Remove adapter contamination and poor-quality bases if present.

```bash
trimmomatic SE -threads 8 -phred33 \
  QC_results/SRR32858437_pass_fastqc.zip \
  FASTQ_files/MCF7_LCOR_OE_trimmed.fastq.gz \
  TRAILING:10
```  
**Tools:** [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic), [`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/), [`fastp`](https://github.com/OpenGene/fastp)

---

## 5. Post-trimming Quality Control
- Re-run FastQC after trimming to ensure cleaning steps were effective.  
**Tool:** [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

---

## 6. Reference Genome Preparation (Indexing)
- Download HISAT2 prebuilt GRCh38 genome index:

```bash

wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz

```
- Download Ensembl GTF annotation:

```bash

wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz

```

---

## 7. Alignment/Mapping
 - Align reads with the Human Genome and convert the SAM file to BAM file.

```bash

hisat2 -q -x grch38/genome -U fastq/SRR7179504_pass.fastq.gz | \
  samtools sort -o alignedreads/LNCAP_Hypoxia_S1.bam
```
**Tools:** [`HISAT2`](https://daehwankimlab.github.io/hisat2/), [`STAR`](https://github.com/alexdobin/STAR), `Bowtie2`

---

## 8. SAM/BAM Processing (Sorting/Indexing)

Convert SAM to BAM; sort and index the alignments.  
**Tools:** [`samtools`](http://www.htslib.org/)

---

## 9. Assessing Alignment Quality

Generate reports on mapping performance and quality.  
**Tools:** `samtools flagstat`, `Qualimap`, `RSeQC`

---

## 10. Feature Counting (Read Quantification)

Count the reads mapping to genes/features.  
**Tools:** [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/), [`HTSeq-count`](https://htseq.readthedocs.io/en/master/)

---

## 11. Downstream Analysis

### 11.1 Data Import, Filtering, and Normalization

- Import feature count matrix into R (or Python)
- Filter out lowly expressed genes
- Normalize counts using methods like TMM (edgeR), DESeq2 normalization (size factors), or TPM/CPM where appropriate

### 11.2 Exploratory Data Analysis (EDA)

- Principal Component Analysis (PCA) to assess sample relationships
- Hierarchical clustering and sample distance heatmaps
- Visualization of library complexity and outliers

### 11.3 Differential Gene Expression Analysis

- Use tools such as [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html), or [`limma-voom`](https://bioconductor.org/packages/release/bioc/html/limma.html)
- Specify design matrix reflecting biological conditions
- Estimate dispersion, fit models, run statistical tests
- Generate lists of differentially expressed genes (DEGs)

### 11.4 Visualization of Results

- MA plots, volcano plots
- Heatmaps for top DEGs across samples
- Gene expression plots (boxplots, barplots)

### 11.5 Functional Enrichment Analysis

- Gene Ontology (GO) enrichment (using `clusterProfiler`, `topGO`, `gProfiler`, etc.)
- Pathway analysis (e.g., KEGG, Reactome)

### 11.6 Advanced/Optional Analyses

- Gene Set Enrichment Analysis (GSEA)
- Splicing analysis (e.g., using rMATS, DEXSeq)
- Cell type deconvolution (if bulk tissue)
- Integration with external data (TCGA, GTEx)
- Visualization in genome browser (IGV)

---

> *Click on any step in [Quick Navigation](#quick-navigation) to jump directly to that section!*
