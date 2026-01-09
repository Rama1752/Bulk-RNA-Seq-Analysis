# ===============================
#    INSTALL AND LOAD PACKAGES
# ===============================

# Install Bioconductor packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "org.Hs.eg.db", "apeglm", "pheatmap", "fgsea", "msigdbr"))
install.packages(c("tidyverse", "ggrepel"))

# Load required libraries
library(DESeq2)
library(apeglm)
library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(msigdbr)

# ===============================
#    DATA IMPORT AND PREPROCESSING
# ===============================

# Read FeatureCounts output
data <- read.table("//wsl.localhost/Ubuntu-24.04/home/rama/RNA-Seq/quants/featurecounts.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

# Remove annotation columns (first 5 columns)
counts <- data[, 6:ncol(data)]

# Rename columns with meaningful sample names
colnames(counts) <- c(
  paste0("Control_rep", 1:5),
  paste0("Cyp_IL1b_rep", 1:5),
  paste0("Cyp_rep", 1:5),
  paste0("IL1b_rep", 1:5)
)

# Filter low-count genes (minimum 10 reads across all samples)
counts_filtered <- counts[rowSums(counts) >= 10, ]

# ===============================
#    CREATE METADATA
# ===============================

conditions <- c(rep("Control", 5),
                rep("Cyp_IL1b", 5),
                rep("Cyp", 5),
                rep("IL1b", 5))

coldata <- data.frame(condition = factor(conditions),
                      row.names = colnames(counts_filtered))

# ===============================
#    BUILD BASE DESEQ2 DATASET
# ===============================

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ condition
)

# ===============================
#    DIFFERENTIAL EXPRESSION ANALYSIS
# ===============================

# ----------------
# Comparison 1: IL1b vs Control
# ----------------

dds1 <- dds
dds1$condition <- relevel(dds1$condition, ref = "Control")
dds1 <- DESeq(dds1)

# Extract results
res1 <- results(dds1, contrast = c("condition", "IL1b", "Control"))

resLFC1 <- lfcShrink(dds1, coef = "condition_IL1b_vs_Control", type = "apeglm")

# Annotate with gene symbols
resLFC1$gene <- mapIds(org.Hs.eg.db,
                       keys = rownames(resLFC1),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Sort by adjusted p-value
resLFC1 <- resLFC1[order(resLFC1$padj, na.last = TRUE), ]

write.csv(as.data.frame(resLFC1), "DESeq/CSV_Files/IL1b_vs_Control_DEGs.csv")

# ----------------
# Comparison 2: Cyp_IL1b vs IL1b
# ----------------

dds2 <- dds
dds2$condition <- relevel(dds2$condition, ref = "IL1b")
dds2 <- DESeq(dds2)

# Extract results
res2 <- results(dds2, contrast = c("condition", "Cyp_IL1b", "IL1b"))
resLFC2 <- lfcShrink(dds2, coef = "condition_Cyp_IL1b_vs_IL1b", type = "apeglm")

# Annotate with gene symbols
resLFC2$gene <- mapIds(org.Hs.eg.db,
                       keys = rownames(resLFC2),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Sort by adjusted p-value
resLFC2 <- resLFC2[order(resLFC2$padj, na.last = TRUE), ]

write.csv(as.data.frame(resLFC2), "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b_DEGs.csv")

# ===============================
#    EXPORT RESULTS FILES
# ===============================

# ----------------
# All genes
# ----------------

# IL1b vs Control
res1_df <- as.data.frame(resLFC1) %>%
  rownames_to_column("ensembl_id")

write.csv(
  res1_df,
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_all_genes.csv",
  row.names = FALSE
)

# Cyp_IL1b vs IL1b
res2_df <- as.data.frame(resLFC2) %>%
  rownames_to_column("ensembl_id")

write.csv(
  res2_df,
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_all_genes.csv",
  row.names = FALSE
)

# ----------------
# Significant genes (padj < 0.05, |log2FC| > 1)
# ----------------

# IL1b vs Control
sig1 <- res1_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(
  sig1,
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_sig_genes.csv",
  row.names = FALSE
)

# Cyp_IL1b vs IL1b
sig2 <- res2_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(
  sig2,
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_sig_genes.csv",
  row.names = FALSE
)

# ----------------
# Top 40 genes (by adjusted p-value)
# ----------------

# IL1b vs Control
top1 <- sig1 %>%
  arrange(padj) %>%
  head(40)

write.csv(
  top1,
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_top40_genes.csv",
  row.names = FALSE
)

# Cyp_IL1b vs IL1b
top2 <- sig2 %>%
  arrange(padj) %>%
  head(40)

write.csv(
  top2,
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_top40_genes.csv",
  row.names = FALSE
)

# ----------------
# Ranked files (for GSEA)
# ----------------

# IL1b vs Control
rank1 <- as.data.frame(res1) %>%
  rownames_to_column("ensembl_id") %>%
  mutate(
    gene = mapIds(
      org.Hs.eg.db,
      keys = ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  ) %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::select(gene, stat) %>%
  dplyr::arrange(desc(stat))

write.table(
  rank1,
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_rank.rnk",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Cyp_IL1b vs IL1b
rank2 <- as.data.frame(res2) %>%
  rownames_to_column("ensembl_id") %>%
  mutate(
    gene = mapIds(
      org.Hs.eg.db,
      keys = ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  ) %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::select(gene, stat) %>%
  dplyr::arrange(desc(stat))

write.table(
  rank2,
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_rank.rnk",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ===============================
#    QUALITY CONTROL VISUALIZATIONS
# ===============================

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# ----------------
# PCA Plot
# ----------------

plot_PCA <- function(vsd.obj) {
  pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"),
         title = "PCA Plot")
}

png("DESeq/Plots/PCA_Plot.png", width = 1200, height = 1000, res = 150)
plot_PCA(vsd)
dev.off()

# ----------------
# Sample Distance Heatmap
# ----------------

plotDists <- function(vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  rownames(sampleDistMatrix) <- colnames(vsd.obj)
  colnames(sampleDistMatrix) <- colnames(vsd.obj)
  
  colors <- colorRampPalette(
    rev(brewer.pal(9, "Blues"))
  )(255)
  
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    main = "Sample-to-Sample Distance Heatmap"
  )
}

png("DESeq/Plots/Distance_Heatmap.png", width = 1200, height = 1000, res = 150)
plotDists(vsd)
dev.off()

# ----------------
# Variable Gene Heatmap
# ----------------

variable_gene_heatmap <- function(vsd.obj,
                                  num_genes = 500,
                                  title = "Top Variable Genes Heatmap") {
  
  # Color palette
  ramp <- colorRampPalette(brewer.pal(11, "RdBu"))
  colors <- ramp(256)[256:1]
  
  # Extract VST counts
  mat <- assay(vsd.obj)
  
  # Compute row variances (genes)
  rv <- rowVars(mat)
  
  # Select top variable genes
  top_mat <- mat[order(rv, decreasing = TRUE)[1:num_genes], ]
  
  # Center genes (row-wise)
  top_mat <- top_mat - rowMeans(top_mat)
  
  # Map Ensembl IDs to gene symbols
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(top_mat),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Replace rownames (keep Ensembl if symbol is NA)
  rownames(top_mat) <- ifelse(
    is.na(gene_symbols),
    rownames(top_mat),
    gene_symbols
  )
  
  # Column annotation (condition only)
  annotation_col <- data.frame(
    condition = colData(vsd.obj)$condition
  )
  rownames(annotation_col) <- colnames(top_mat)
  
  # Plot heatmap
  pheatmap(
    top_mat,
    color = colors,
    annotation_col = annotation_col,
    fontsize_col = 8,
    fontsize_row = 250 / num_genes,
    border_color = NA,
    main = title
  )
}

png("DESeq/Plots/variable_gene_Heatmap.png", width = 1200, height = 1000, res = 150)
variable_gene_heatmap(vsd, num_genes = 40,
                      title = "Top 40 Variable Genes")
dev.off()

# ===============================
#    DIFFERENTIAL EXPRESSION VISUALIZATIONS
# ===============================

# ----------------
# Volcano Plots
# ----------------

volcano_plot_gg <- function(csv_file, title_text, label_n = 10) {
  
  # Read CSV
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Basic filtering and transformation
  df <- df %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = ifelse(padj < 0.05 & abs(log2FoldChange) >= 1,
                   "Significant", "Not significant"),
      neg_log10_padj = -log10(padj)
    )
  
  # Genes to label (top by padj)
  label_df <- df %>%
    filter(sig == "Significant") %>%
    arrange(padj) %>%
    head(label_n)
  
  ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Not significant" = "blue",
                                  "Significant" = "red")) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", linewidth = 0.5) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = title_text,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# IL1b vs Control
png("DESeq/Plots/Volcano_plot_IL1b_vs_Control.png", width = 1200, height = 1000, res = 150)
volcano_plot_gg(
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_all_genes.csv",
  "Volcano Plot: IL1b vs Control"
)
dev.off()

# Cyp_IL1b vs IL1b
png("DESeq/Plots/Volcano_plot_Cyp_IL1b_vs_IL1b.png", width = 1200, height = 1000, res = 150)
volcano_plot_gg(
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_all_genes.csv",
  "Volcano Plot: Cyp_IL1b vs IL1b"
)
dev.off()

# ----------------
# Log2FoldChange Comparison Plot
# ----------------

# Read results files
res1 <- read.csv(
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_all_genes.csv",
  header = TRUE
)

res2 <- read.csv(
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_all_genes.csv",
  header = TRUE
)

compare_significant_genes <- function(res1, res2, padj_cutoff = 0.0001, 
                                      ngenes = 250, nlabel = 10, 
                                      samplenames = c("comparison1", "comparison2"), 
                                      title = "") {
  # Get list of most upregulated or downregulated genes for each results table
  genes1 <- rbind(head(res1[which(res1$padj < padj_cutoff), ], ngenes), 
                  tail(res1[which(res1$padj < padj_cutoff), ], ngenes))
  genes2 <- rbind(head(res2[which(res2$padj < padj_cutoff), ], ngenes), 
                  tail(res2[which(res2$padj < padj_cutoff), ], ngenes))
  
  # Combine the data from both tables
  de_union <- union(genes1$ensembl_id, genes2$ensembl_id)
  res1_union <- res1[match(de_union, res1$ensembl_id), ][c("ensembl_id", "log2FoldChange", "gene")]
  res2_union <- res2[match(de_union, res2$ensembl_id), ][c("ensembl_id", "log2FoldChange", "gene")]
  combined <- left_join(res1_union, res2_union, by = "ensembl_id", suffix = samplenames)
  
  # Identify overlap between genes in both tables
  combined$de_condition <- NA_character_
  combined$de_condition[which(combined$ensembl_id %in% intersect(genes1$ensembl_id, genes2$ensembl_id))] <- "Significant in Both"
  combined$de_condition[which(combined$ensembl_id %in% setdiff(genes1$ensembl_id, genes2$ensembl_id))] <- paste0("Significant in ", samplenames[1])
  combined$de_condition[which(combined$ensembl_id %in% setdiff(genes2$ensembl_id, genes1$ensembl_id))] <- paste0("Significant in ", samplenames[2])
  
  # Find the top most genes within each condition to label on the graph
  label1 <- rbind(head(combined[which(combined$de_condition == paste0("Significant in ", samplenames[1])), ], nlabel),
                  tail(combined[which(combined$de_condition == paste0("Significant in ", samplenames[1])), ], nlabel))
  label2 <- rbind(head(combined[which(combined$de_condition == paste0("Significant in ", samplenames[2])), ], nlabel),
                  tail(combined[which(combined$de_condition == paste0("Significant in ", samplenames[2])), ], nlabel))
  label3 <- rbind(head(combined[which(combined$de_condition == "Significant in Both"), ], nlabel),
                  tail(combined[which(combined$de_condition == "Significant in Both"), ], nlabel))
  combined_labels <- rbind(label1, label2, label3)
  
  # Plot the genes based on log2FoldChange, color coded by significance
  ggplot(
    combined,
    aes_string(
      x = paste0("`log2FoldChange", samplenames[1], "`"),
      y = paste0("`log2FoldChange", samplenames[2], "`")
    )
  ) +
    geom_point(aes(color = de_condition), size = 0.7) +
    scale_color_manual(values = c("#00BA38", "#619CFF", "#F8766D")) +
    geom_text_repel(
      data = combined_labels,
      aes_string(
        label = paste0("`gene", samplenames[1], "`"),
        color = "de_condition"
      ),
      show.legend = FALSE,
      size = 3
    ) +
    geom_vline(xintercept = 0, size = 0.3, linetype = 2) +
    geom_hline(yintercept = 0, size = 0.3, linetype = 2) +
    labs(
      title = title,
      x = paste0("log2FoldChange in ", samplenames[1]),
      y = paste0("log2FoldChange in ", samplenames[2])
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

png("DESeq/Plots/Log2FoldChange_Comparison_Plot.png", width = 1200, height = 1000, res = 150)
compare_significant_genes(
  res1,
  res2,
  samplenames = c("IL1b_vs_Control", "Cyp_IL1b_vs_IL1b"),
  title = "Gene-level reversal of IL1b effects by Cyp"
)
dev.off()

# ----------------
# DE Gene Heatmaps 
# ----------------

DE_gene_heatmap <- function(
    res,
    vsd,
    padj_cutoff = 0.05,
    ngenes = 30,
    title = "Top DE genes",
    condition_order,
    cluster_cols = FALSE
) {
  
  # Select significant genes
  sig_genes <- res %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::arrange(desc(abs(log2FoldChange))) %>%
    head(ngenes)
  
  gene_ids <- rownames(sig_genes)
  
  # Extract VST expression
  mat <- assay(vsd)[gene_ids, ]
  
  # Replace Ensembl IDs with symbols
  rownames(mat) <- ifelse(
    is.na(sig_genes$gene),
    gene_ids,
    sig_genes$gene
  )
  
  # Sample annotation
  annotation_col <- data.frame(
    Condition = colData(vsd)$condition
  )
  rownames(annotation_col) <- colnames(mat)
  
  # Enforce comparison-specific order (KEY FIX)
  ord <- order(factor(annotation_col$Condition, levels = condition_order))
  mat <- mat[, ord]
  annotation_col <- annotation_col[ord, , drop = FALSE]
  
  # Color palette
  colors <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(9, "RdBu"))
  )(255)
  
  # Plot heatmap
  pheatmap::pheatmap(
    mat,
    color = colors,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = cluster_cols,
    annotation_col = annotation_col,
    fontsize_row = 200 / ngenes,
    fontsize_col = 9,
    border_color = NA,
    main = title
  )
}

png("DESeq/Plots/DE_gene_heatmap_IL1b_vs_Control.png", width = 1200, height = 1000, res = 150)
DE_gene_heatmap(
  res = resLFC1,
  vsd = vsd,
  padj_cutoff = 0.05,
  ngenes = 30,
  title = "Top IL-1β-responsive genes",
  condition_order = c("Cyp", "Control", "IL1b", "Cyp_IL1b"),
  cluster_cols = FALSE
)
dev.off()

png("DESeq/Plots/DE_gene_heatmap_Cyp_IL1b_vs_IL1b.png", width = 1200, height = 1000, res = 150)
DE_gene_heatmap(
  res = resLFC2,
  vsd = vsd,
  padj_cutoff = 0.05,
  ngenes = 30,
  title = "Suppression of IL-1β–induced genes by Cyp",
  condition_order = c("Control", "IL1b", "Cyp_IL1b", "Cyp"),
  cluster_cols = FALSE
)
dev.off()

# ===============================
#    GENE SET ENRICHMENT ANALYSIS (GSEA)
# ===============================

# ----------------
# Load Hallmark Gene Sets
# ----------------

hallmark_sets <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
)

hallmark_list <- split(
  hallmark_sets$gene_symbol,
  hallmark_sets$gs_name
)

# ----------------
# Prepare Ranked Gene Lists
# ----------------

# IL1b vs Control
rank1 <- read.table(
  "DESeq/CSV_Files/IL1b_vs_Control/IL1b_vs_Control_rank.rnk",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

gene_list_IL1b <- rank1$stat
names(gene_list_IL1b) <- rank1$gene

# Remove NA and duplicated gene symbols
gene_list_IL1b <- gene_list_IL1b[!is.na(names(gene_list_IL1b))]
gene_list_IL1b <- gene_list_IL1b[!duplicated(names(gene_list_IL1b))]

# Ensure sorted
gene_list_IL1b <- sort(gene_list_IL1b, decreasing = TRUE)

# Cyp_IL1b vs IL1b
rank2 <- read.table(
  "DESeq/CSV_Files/Cyp_IL1b_vs_IL1b/Cyp_IL1b_vs_IL1b_rank.rnk",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


gene_list_Cyp <- rank2$stat
names(gene_list_Cyp) <- rank2$gene

gene_list_Cyp <- gene_list_Cyp[!is.na(names(gene_list_Cyp))]

gene_list_Cyp <- gene_list_Cyp[!duplicated(names(gene_list_Cyp))]
gene_list_Cyp <- sort(gene_list_Cyp, decreasing = TRUE)

# ----------------
# Run fgsea
# ----------------

fgsea_IL1b <- fgsea(
  pathways = hallmark_list,
  stats    = gene_list_IL1b,
  minSize  = 15,
  maxSize  = 500,
  nperm    = 10000
)

fgsea_Cyp <- fgsea(
  pathways = hallmark_list,
  stats    = gene_list_Cyp,
  minSize  = 15,
  maxSize  = 500,
  nperm    = 10000
)

# ----------------
# Compare GSEA Results
# ----------------

compare_hallmark <- fgsea_IL1b %>%
  dplyr::select(pathway, NES_IL1b = NES, padj_IL1b = padj) %>%
  dplyr::left_join(
    fgsea_Cyp %>%
      dplyr::select(pathway, NES_Cyp = NES, padj_Cyp = padj),
    by = "pathway"
  )

# View specific pathways of interest
compare_hallmark %>%
  filter(grepl("INFLAMMATORY|IL6|TNFA|CYTOKINE|AUTOPHAGY", pathway))

# ----------------
# Enrichment Plot
# ----------------

png("DESeq/Plots/HALLMARK_INFLAMMATORY_RESPONSE_enrichment_IL1b_vs_Control.png", width = 1200, height = 1000, res = 150)
plotEnrichment(
  hallmark_list[["HALLMARK_INFLAMMATORY_RESPONSE"]],
  gene_list_IL1b
) +
  labs(title = "Inflammatory Response – IL1b vs Control")
dev.off()

png("DESeq/Plots/HALLMARK_INFLAMMATORY_RESPONSE_enrichment_Cyp_IL1b_vs_IL1b.png", width = 1200, height = 1000, res = 150)
plotEnrichment(
  hallmark_list[["HALLMARK_INFLAMMATORY_RESPONSE"]],
  gene_list_Cyp
) +
  labs(title = "Inflammatory Response – Cyp_IL1b vs IL1b")
dev.off()

# ----------------
# Waterfall Plots
# ----------------

waterfall_plot <- function(fgsea_results, graph_title) {
  
  fgsea_results %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::mutate(
      short_name = stringr::str_replace(pathway, "HALLMARK_", ""),
      regulation = ifelse(NES > 0, "Activated", "Suppressed")
    ) %>%
    ggplot(aes(x = reorder(short_name, NES), y = NES)) +
    geom_col(aes(fill = regulation)) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Activated"  = "#D55E00",  # orange
        "Suppressed" = "#0072B2"   # blue
      )
    ) +
    labs(
      title = graph_title,
      x = "Hallmark Pathway",
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.title = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
    )
}


png("DESeq/Plots/Waterfall_plot_IL1b_vs_Control.png", width = 1200, height = 1000, res = 150)
waterfall_plot(
  fgsea_IL1b,
  "Hallmark pathways altered by IL1b treatment"
)
dev.off()

png("DESeq/Plots/Waterfall_plot_Cyp_IL1b_vs_IL1b.png", width = 1200, height = 1000, res = 150)
waterfall_plot(
  fgsea_Cyp,
  "Hallmark pathways altered by Cyproheptadine under IL1b"
)
dev.off()


