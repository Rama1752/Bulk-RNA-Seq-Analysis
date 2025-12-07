# Bulk RNA-Seq Analysis

🎯 Project Overview:

This project demonstrates a complete Bulk RNA-Seq workflow analyzing how human cartilage cells (chondrocytes) respond to cyproheptadine, a drug that shows promise for treating osteoarthritis. We studied gene expression changes to understand how this drug affects cartilage cells at the molecular level.

The study examined the effects of:
- **Cyproheptadine** - a drug candidate that may help protect cartilage and slow down osteoarthritis
- **IL-1β** - an inflammatory molecule that damages cartilage (mimics what happens in osteoarthritis)
- **Combination treatment** - cyproheptadine given together with IL-1β to see if the drug can protect cells from inflammation

This dataset helps us compare how cyproheptadine works in healthy cells versus inflamed cells, giving us insights into how it might treat osteoarthritis by protecting the cartilage in our joints.

## 🧪 Sample Information

| Sample ID      | Treatment Condition | SRA Accession                                                      |
|----------------|---------------------|--------------------------------------------------------------------|
| Sample 1 - 5   | Cyp + IL-1β         | SRR32684363, SRR32684364, SRR32684365, SRR32684366, SRR32684367    |
| Sample 6 - 10  | IL-1β               | SRR32684368, SRR32684369, SRR32684370, SRR32684371, SRR32684372    |
| Sample 11 - 15 | Cyproheptadine      | SRR32684373, SRR32684374, SRR32684375, SRR32684376, SRR32684377    |
| Sample 16 - 20 | DMSO (Control)      | SRR32684378, SRR32684379, SRR32684380, SRR32684381, SRR32684382    |

📂 Dataset Information:
- **GEO Accession:** [GSE291878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291878)
- **Experiment Type:** Expression profiling by high-throughput sequencing
- **Description:** RNA-Seq experiment in human cartilage cells to see how cyproheptadine changes gene activity in normal conditions and during inflammation. Osteoarthritis is a common joint disease with no cure, and this study helps identify cyproheptadine as a potential treatment option.
- **Citations:**
  1. Kurakazu I, Olmer M, Swahn H, Myers K et al. (2025) Histamine H1 receptor inverse agonists improve structure and pain in an osteoarthritis mouse model. J Clin Invest 135(21). [PMID: 40875519]

## Quick Navigation

- [1. Downloading necessary tools and make directories](#1-downloading-necessary-tools-and-make-directories)
- [2. Fetching SRA Files](#2-fetching-sra-files)
- [3. FastQC - Quality Control](#3-fastqc---quality-control)
- [4. MultiQC](#4-multiqc)
- [5. Trimming Reads (optional)](#5-trimming-reads-optional)
- [6. Post-trimming Quality Control](#6-post-trimming-quality-control)
- [7. Reference Genome Preparation (Indexing)](#7-reference-genome-preparation-indexing)
- [8. Alignment/Mapping](#8-alignmentmapping)
- [9. BAM index file](#9-bam-index-file)
- [10. Assessing Alignment Quality](#10-assessing-alignment-quality)
- [11. Convert GTF to BED](#11-convert-gtf-to-bed)
- [12. Determine Library Strandedness](#12-determine-library-strandedness)
- [13. Feature Counting (Read Quantification)](#13-feature-counting-read-quantification)
- [14. Downstream Analysis](#14-downstream-analysis)

---

## Workflow

---

## 1. Downloading necessary tools and make directories
- Download all the necessary tools and make directories.

```bash

conda install -y -c bioconda -c conda-forge fastqc multiqc sra-tools hisat2 samtools cutadapt subread qualimap rseqc bedops

mkdir -p SRA_files FASTQ_files FASTQC_reports Multiqc_reports reference aligned_reads quants rnaseq_qc_results trimmed trimmed_qc_reports

```
---

## 2. Fetching SRA Files
- Download SRA files using the SRA toolkit.
- Convert the SRA files into FASTQ files using fastq-dump.
- Save the file in .gz format i.e. gzip file so that it would take less space.

```bash
#Download SRA files
prefetch SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512 SRR1039513 \
SRR1039514 SRR1039515 SRR1039516 SRR1039517 SRR1039518 SRR1039519 SRR1039520 \
SRR1039521 SRR1039522 SRR1039523 --progress

#Covert SRA files to FASTQ files
fastq-dump --outdir FASTQ_files --gzip --skip-technical \
--readids --read-filter pass --dumpbase --split-3 --clip \
SRR1039508/SRR1039508.sra
```
---

## 3. FastQC - Quality Control
- Check the quality of raw sequencing reads using FastQC.
- Add threads as per the number of CPU cores available.
```bash

#Run FASTQC
fastqc FASTQ_files/*.fastq.gz -o FASTQC_reports/ --threads 8

```
---

## 4. MultiQC 
- Merge all the fastqc reports to get a summarised report of it.

```bash

multiqc FASTQC_reports/ -o Multiqc_reports

```

## 5. Trimming Reads
- Remove adapter contamination and poor-quality bases if present.
- Copy all the over-represented seq from fastqc reports for trimming.
- -a seq is of R1 read and -A is of R2 read.
- Add threads -j as per the number of CPU cores available.

```bash

cutadapt \
  -j 8 \
  -a ACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTT \
  -a ACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTT \
  -a ACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTT \
  -a CACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTG \
  -A GTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA \
  -A GTCGTGTAGGGAAAGAGGGTAGATCTCGGTGGTCGCCGTATCATTAAAAA \
  -A CGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAA \
  --minimum-length 36 \
  -o trimmed/SRR1039508_pass_1P.cutadapt.fastq.gz \
  -p trimmed/SRR1039508_pass_2P.cutadapt.fastq.gz \
  FASTQ_files/SRR1039508_pass_1.fastq.gz \
  FASTQ_files/SRR1039508_pass_2.fastq.gz \
  > trimmed/SRR1039508_cutadapt_report.txt 2>&1

```  
---

## 6. Post-trimming Quality Control
- Re-run FastQC after trimming to ensure cleaning steps were effective.

```bash

fastqc trimmed/* -o trimmed_qc_reports

```

---

## 7. Reference Genome Preparation (Indexing)
- Download HISAT2 prebuilt GRCh37 genome index and ensemble gtf annotation.

```bash

wget https://genome-idx.s3.amazonaws.com/hisat/grch37_genome.tar.gz
tar -xvzf grch37_genome.tar.gz -C reference/

```
- Download Ensembl GTF annotation:

```bash

wget -P reference/ https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz

```

---

## 8. Alignment/Mapping
 - Rename the files for better understanding.
 - Align reads with the Human Genome and convert the SAM file to BAM file.
 - Add threads -p or -@ as per the number of CPU cores available.

```bash

#Renaming the files:
mv trimmed/SRR1039508_pass_1P.cutadapt.fastq.gz trimmed/N61311_untreated_R1.fastq.gz
mv trimmed/SRR1039508_pass_2P.cutadapt.fastq.gz trimmed/N61311_untreated_R2.fastq.gz
mv trimmed/SRR1039509_pass_1P.cutadapt.fastq.gz trimmed/N61311_Dex_R1.fastq.gz
mv trimmed/SRR1039509_pass_2P.cutadapt.fastq.gz trimmed/N61311_Dex_R2.fastq.gz
.....

#Aligning the fastq files with genome
hisat2 -p 6 -q -x reference/grch37/genome -1 trimmed/N61311_untreated_R1.fastq.gz -2 trimmed/N61311_untreated_R2.fastq.gz | \
  samtools sort -@ 4 -o aligned_reads/N61311_untreated.bam

```
---

## 9. BAM index file
- Creates .bai index file for fast random access.

```bash

samtools index -@ 8 -M aligned_reads/*.bam

```
---

## 10. Assessing Alignment Quality
- Generate reports on mapping performance and quality.

```bash

#QC Check
qualimap rnaseq -bam aligned_reads/MDA_MB_231_LCOR_OE.bam -gtf reference/Homo_sapiens.GRCh38.115.gtf \
 -outdir rnaseq_qc_results/MDA_MB_231_LCOR_OE --java-mem-size=10G
 
```

---

## 11. Convert GTF to BED
- Create a .bed file from Human .gtf file which is required to check strandedness of RNA-Seq data using RSeQC.

```bash

gtf2bed < reference/Homo_sapiens.GRCh38.115.gtf > reference/Homo_sapiens.GRCh38.115.bed

```
---

## 12. Determine Library Strandedness
- RSeQC is used to check the strandedness of the RNA-Seq data.
- Checking strandedness is necessary cause it assigns the reads to the correct gene, when genes overlap on opposite strands.

```bash

infer_experiment.py -i aligned_reads/MDA_MB_231_LCOR_OE.bam \
  -r reference/Homo_sapiens.GRCh38.115.bed

```

---
## 13. Feature Counting (Read Quantification)
- Count the reads mapping to genes/features.  

```bash

featureCounts -S 2 -a reference/Homo_sapiens.GRCh38.115.gtf \
  -o quants/featurecounts.txt aligned_reads/*.bam

```

---

## 14. Downstream Analysis

### 14.1 Data Import, Filtering, and Normalization

- Import feature count matrix into R (or Python)
- Filter out lowly expressed genes
- Normalize counts using methods like TMM (edgeR), DESeq2 normalization (size factors), or TPM/CPM where appropriate

### 14.2 Exploratory Data Analysis (EDA)

- Principal Component Analysis (PCA) to assess sample relationships
- Hierarchical clustering and sample distance heatmaps
- Visualization of library complexity and outliers

### 14.3 Differential Gene Expression Analysis

- Use tools such as [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html), or [`limma-voom`](https://bioconductor.org/packages/release/bioc/html/limma.html)
- Specify design matrix reflecting biological conditions
- Estimate dispersion, fit models, run statistical tests
- Generate lists of differentially expressed genes (DEGs)

### 14.4 Visualization of Results

- MA plots, volcano plots
- Heatmaps for top DEGs across samples
- Gene expression plots (boxplots, barplots)

### 14.5 Functional Enrichment Analysis

- Gene Ontology (GO) enrichment (using `clusterProfiler`, `topGO`, `gProfiler`, etc.)
- Pathway analysis (e.g., KEGG, Reactome)

### 14.6 Advanced/Optional Analyses

- Gene Set Enrichment Analysis (GSEA)
- Splicing analysis (e.g., using rMATS, DEXSeq)
- Cell type deconvolution (if bulk tissue)
- Integration with external data (TCGA, GTEx)
- Visualization in genome browser (IGV)

---

> *Click on any step in [Quick Navigation](#quick-navigation) to jump directly to that section!*
