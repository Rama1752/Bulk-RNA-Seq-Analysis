# RNA-Seq Analysis Pipeline Outline

## Quick Navigation

- [1. Fetching SRA Files](#1-fetching-sra-files)
- [2. FastQC - Quality Control](#2-fastqc-quality-control)
- [3. Trimming Reads](#3-trimming-reads)
- [4. Post-trimming Quality Control](#4-post-trimming-quality-control)
- [5. Reference Genome Preparation (Indexing)](#5-reference-genome-preparation-indexing)
- [6. Alignment/Mapping](#6-alignment-mapping)
- [7. SAM/BAM Processing (Sorting/Indexing)](#7-sam-bam-processing-sorting-indexing)
- [8. Assessing Alignment Quality](#8-assessing-alignment-quality)
- [9. Feature Counting (Read Quantification)](#9-feature-counting-read-quantification)
- [10. Downstream Analysis](#10-downstream-analysis)

---

## 1. Fetching SRA Files

Download sample FASTQ files from the Sequence Read Archive (SRA).  
**Tools:** `sra-tools` (`prefetch`, `fasterq-dump`)

---

## 2. FastQC - Quality Control

Assess the quality of raw sequencing reads using FastQC.  
**Tool:** [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

---

## 3. Trimming Reads

Remove adapter contamination and poor-quality bases.  
**Tools:** [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic), [`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/), [`fastp`](https://github.com/OpenGene/fastp)

---

## 4. Post-trimming Quality Control

Re-run QC after trimming to ensure cleaning steps were effective.  
**Tool:** `FastQC`

---

## 5. Reference Genome Preparation (Indexing)

Index the reference genome to allow quick alignment.  
**Tools:** `HISAT2-build`, `STAR`, `Bowtie2-build`, etc.

---

## 6. Alignment/Mapping

Map the reads to the reference genome/transcriptome.  
**Tools:** [`HISAT2`](https://daehwankimlab.github.io/hisat2/), [`STAR`](https://github.com/alexdobin/STAR), `Bowtie2`

---

## 7. SAM/BAM Processing (Sorting/Indexing)

Convert SAM to BAM; sort and index the alignments.  
**Tools:** [`samtools`](http://www.htslib.org/)

---

## 8. Assessing Alignment Quality

Generate reports on mapping performance and quality.  
**Tools:** `samtools flagstat`, `Qualimap`, `RSeQC`

---

## 9. Feature Counting (Read Quantification)

Count the reads mapping to genes/features.  
**Tools:** [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/), [`HTSeq-count`](https://htseq.readthedocs.io/en/master/)

---

## 10. Downstream Analysis

### 10.1 Data Import, Filtering, and Normalization

- Import feature count matrix into R (or Python)
- Filter out lowly expressed genes
- Normalize counts using methods like TMM (edgeR), DESeq2 normalization (size factors), or TPM/CPM where appropriate

### 10.2 Exploratory Data Analysis (EDA)

- Principal Component Analysis (PCA) to assess sample relationships
- Hierarchical clustering and sample distance heatmaps
- Visualization of library complexity and outliers

### 10.3 Differential Gene Expression Analysis

- Use tools such as [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html), or [`limma-voom`](https://bioconductor.org/packages/release/bioc/html/limma.html)
- Specify design matrix reflecting biological conditions
- Estimate dispersion, fit models, run statistical tests
- Generate lists of differentially expressed genes (DEGs)

### 10.4 Visualization of Results

- MA plots, volcano plots
- Heatmaps for top DEGs across samples
- Gene expression plots (boxplots, barplots)

### 10.5 Functional Enrichment Analysis

- Gene Ontology (GO) enrichment (using `clusterProfiler`, `topGO`, `gProfiler`, etc.)
- Pathway analysis (e.g., KEGG, Reactome)

### 10.6 Advanced/Optional Analyses

- Gene Set Enrichment Analysis (GSEA)
- Splicing analysis (e.g., using rMATS, DEXSeq)
- Cell type deconvolution (if bulk tissue)
- Integration with external data (TCGA, GTEx)
- Visualization in genome browser (IGV)

---

> *Click on any step in [Quick Navigation](#quick-navigation) to jump directly to that section!*
