# Bulk RNA-Seq Analysis

A comprehensive pipeline for analyzing bulk RNA-sequencing data with quality control, alignment, quantification, and differential expression analysis.

---

## Table of Contents

- [Overview](#overview)
- [Dataset Information](#dataset-information)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Key Findings Summary](#key-findings-summary)
- [Outputs](#outputs)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

This project provides an end-to-end pipeline for bulk RNA-Seq analysis, including:
- Quality control and preprocessing
- Read alignment and quantification
- Differential expression analysis
- Visualization and interpretation

---

## Dataset Information

### Source Details
- **Organism**: *Homo sapiens* (Human)
- **Tissue Type**: Various tissues
- **Sample Count**: Multiple biological replicates
- **Sequencing Platform**: Illumina sequencing
- **Read Type**: Paired-end reads

### Data Specifications
- **Quality Control**: FastQC analysis performed
- **Alignment Reference**: Human genome (GRCh38)
- **Quantification Method**: Transcript-level abundance estimation

---

## Features

- ✅ **Quality Control**: FastQC reports and quality metrics
- ✅ **Read Alignment**: STAR alignment to reference genome
- ✅ **Expression Quantification**: Kallisto/Salmon for transcript abundance
- ✅ **Differential Expression**: DESeq2 analysis
- ✅ **Visualization**: PCA, heatmaps, volcano plots, and MA plots
- ✅ **Statistical Analysis**: Comprehensive statistical testing
- ✅ **Gene Annotation**: Functional annotation of results

---

## Installation

### Prerequisites
- R (≥ 4.0)
- Python (≥ 3.7)
- STAR aligner
- Kallisto or Salmon

### Setup Steps

```bash
# Clone the repository
git clone https://github.com/Rama1752/Bulk-RNA-Seq-Analysis.git
cd Bulk-RNA-Seq-Analysis

# Install R dependencies
Rscript install_dependencies.R

# Install Python packages
pip install -r requirements.txt
```

---

## Usage

### Running the Complete Pipeline

```bash
# Execute the main analysis script
bash run_analysis.sh --input data/ --output results/

# Or use the R script directly
Rscript main_analysis.R --input data/ --output results/
```

### Individual Steps

#### Step 1: Quality Control
```bash
fastqc data/*.fastq.gz -o qc_reports/
```

#### Step 2: Alignment
```bash
STAR --genomeDir reference/ \
     --readFilesIn data/sample_R1.fastq.gz data/sample_R2.fastq.gz \
     --outFileNamePrefix results/sample_
```

#### Step 3: Quantification
```bash
kallisto quant -i reference/transcripts.idx \
               -o results/sample_quant \
               data/sample_R1.fastq.gz data/sample_R2.fastq.gz
```

#### Step 4: Differential Expression Analysis
```bash
Rscript differential_expression.R --input results/ --output de_results/
```

---

## Results

### Analysis Overview

The pipeline generates comprehensive results across multiple analysis stages, from quality assessment through statistical testing and visualization.

### Key Findings Summary

#### Expression Patterns
- Identified significant expression variation across samples
- Detected tissue-specific gene expression signatures
- Observed consistent clustering patterns in dimensionality reduction

#### Differential Expression Highlights
- **Total Genes Tested**: Comprehensive genome-wide analysis
- **Significant Genes**: Multiple genes with adjusted p-value < 0.05
- **Effect Sizes**: Substantial log2 fold-changes detected in key genes
- **Biological Significance**: Results validated through functional enrichment

#### Quality Metrics
- High-quality sequencing across all samples
- Consistent alignment rates across replicates
- Low technical variation within biological groups

### Data Visualization

#### Overview Plots
- **PCA Analysis**: Principal component analysis revealing sample clustering
- **Heatmaps**: Expression patterns across genes and samples
- **Distribution Plots**: Library size and quality score distributions

#### Differential Expression Visualizations
- **Volcano Plots**: Statistical significance vs. effect size
- **MA Plots**: Mean abundance vs. log fold-change
- **Expression Profiles**: Individual gene expression patterns

---

## Outputs

### Directory Structure

```
results/
├── qc/
│   ├── fastqc_reports/
│   │   └── *.html
│   └── quality_metrics.csv
│
├── alignment/
│   ├── bam_files/
│   │   └── *.bam
│   └── alignment_stats.txt
│
├── quantification/
│   ├── abundance_estimates/
│   │   └── abundance.h5
│   └── tpm_matrix.csv
│
├── differential_expression/
│   ├── de_results.csv
│   ├── gene_annotations.csv
│   └── enrichment_analysis.xlsx
│
└── visualizations/
    ├── pca_plot.pdf
    ├── heatmap.pdf
    ├── volcano_plot.pdf
    ├── ma_plot.pdf
    └── expression_profiles/
```

### Output Files

#### Quality Control Files
- `quality_metrics.csv`: Summary statistics for all samples
- `fastqc_reports/`: Detailed quality reports in HTML format

#### Expression Data
- `tpm_matrix.csv`: Transcript-level quantification
- `abundance_estimates/`: Kallisto/Salmon abundance data

#### Differential Expression Results
- `de_results.csv`: Complete statistical testing results with:
  - Gene identifiers
  - Log2 fold-changes
  - P-values and adjusted p-values
  - Expression levels
  - Significance flags

#### Visualizations
- `pca_plot.pdf`: Sample clustering analysis
- `heatmap.pdf`: Gene expression heatmaps
- `volcano_plot.pdf`: Differential expression volcano plot
- `ma_plot.pdf`: Bland-Altman style MA plot
- `expression_profiles/`: Individual gene expression patterns

#### Functional Analysis
- `enrichment_analysis.xlsx`: Gene ontology and pathway enrichment results
- `gene_annotations.csv`: Functional annotations for significant genes

---

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

### Development Guidelines

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m 'Add your feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

**Last Updated**: 2026-01-08  
**Maintainer**: Rama1752
