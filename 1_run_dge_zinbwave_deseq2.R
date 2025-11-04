#!/usr/bin/env Rscript

library(BiocManager)
library(zinbwave)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(tidyverse)
library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)
library(tximport)
library(tximportData)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(dplyr)
library(scRNAseq)
library(sparseMatrixStats)
library(edgeR)
library(limma)
library(DESeq2)


# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())

# INPUT ----

# Input 
projectdirectory  <- "./RNASeq_Analysis"
input_dir         <- file.path(projectdirectory,"input")
resource_dir      <- file.path(projectdirectory,"resources")
processed_gene_counts_file    <- file.path(input_dir,"input_gene_count_matrix.csv")
metadata_tbm_file             <- file.path(input_dir, "metadata_TBM_formatted.csv")
biomart_metadata_fh           <- file.path(resource_dir,"biomart_protein_coding_genes.tsv")

# Output DGE 
output_dir            <- file.path(projectdirectory,"tables","DGE_results")
dir.create(output_dir)

# DESeq2 output
# - ZI
dge_results_file1 <- file.path(output_dir,"dge_results_zi_deseq.csv")
dge_results_file1_sig <- file.path(output_dir,"dge_results_zi_deseq_significant.csv")
bcv_results_fh        <- file.path(output_dir,"BCV_table_zi_deseq.csv")
# non-survivor vs survivor
dge_results_file1_rev <- file.path(output_dir,"dge_results_zi_deseq_NonSurvivor.csv")
dge_results_file1_sig_rev <- file.path(output_dir,"dge_results_zi_deseq_significant_NonSurvivor.csv")


# 1. Load dataframes ----

### a. load input dataframes ----

# ensembl gene counts
df <- read.csv(processed_gene_counts_file, row.names = 1,check.names = F)

# TBM metadata
metadata            <- read.csv(metadata_tbm_file, stringsAsFactors = FALSE,check.names = F)
row.names(metadata) <- metadata$Sample_ID
metadata            <- metadata[metadata$Group_rx %in% "1",]

# Overlap samples
metadata <- metadata[metadata$Sample_ID %in% names(df),]
metadata$Mortality <- as.factor(metadata$short_term_Mortality)
metadata$Gender    <- as.factor(metadata$gender)
metadata$Batch     <- as.numeric(metadata$BatchAnalysis)

# Calculate Non-Zero genes
metadata$Cohort <- "TBM"
metadata$totalcounts <- colSums(df[,metadata$Sample_ID])
metadata$numNONZero <- as.numeric(colSums(df[,metadata$Sample_ID] != 0, na.rm = T)) # number of non-zero genes
metadata$cohortbatch <- paste(metadata$Cohort, metadata$Batch, sep = "_")
metadata$mortality <-  ifelse(as.numeric(metadata$Mortality) > 1, "NonSurvivor", "Survivor")

# 2. Filter out genes ----
# * Minimum QC is to require non-zero counts in at least 20% of samples per gene
CUTOFF = 0.2*ncol(df)

# Filter for genes detected in >=20% of samples
keep <- rowSums(df>0) >= CUTOFF
keep <- names(keep[keep == TRUE])
total_genes_kept <- length(keep) # number of genes retained

# Filter DGE object
df_filt <- df[keep,metadata$Sample_ID]
df_filt <- na.omit(df_filt)


# 3. Filter for protein coding genes
pc_genes        <- read.delim(biomart_metadata_fh, stringsAsFactors = FALSE,check.names = F)
df_filt$ensembl <- tstrsplit(row.names(df_filt),"_")[[2]]
df_fin          <- subset(df_filt, ensembl %in% pc_genes$ensembl)



# Keep gene lengths and GC content
gene_meta            <- pc_genes[pc_genes$ensembl %in% df_fin$ensembl,]
row.names(gene_meta) <- gene_meta$ensembl
gene_meta            <- gene_meta[df_fin$ensembl,]
row.names(gene_meta) <- row.names(df_fin)
gene_meta            <- gene_meta[,c("length","gccontent")]


df_filt$ensembl <- NULL
df_fin$ensembl  <- NULL

# Convert to matrix
df_fin <- as.matrix(df_fin)

# 3. Make object ----
# * Remove rows with all zeros
se_obj <- SummarizedExperiment(
  assays = list(counts = df_fin), colData = data.frame(metadata), rowData = gene_meta
)

#filter out the rows that are all zeros (should be no rows filtered due to 20% threshold)
filter <- rowSums(assay(se_obj))
se_obj <- se_obj[filter != 0, ]


# 4. Add Quartile categories ----


# 5. Modeling and  calculate weights ----
cat(paste0(">>> Zinbwave dimension reduction \n\n"))

MODEL_COL <- model.matrix(~cohortbatch, data=colData(se_obj) )
zinb_dim_reduc <- zinbwave(se_obj, K = 2, X=MODEL_COL, BPPARAM=BiocParallel::SerialParam(), zeroinflation = TRUE, verbose=FALSE, observationalWeights = TRUE )
weights        <- assay(zinb_dim_reduc, "weights")


# 6. Setup DGE input object and model ----
counts(zinb_dim_reduc)        <- as.matrix(counts(zinb_dim_reduc))
dds <- DESeqDataSet(zinb_dim_reduc, design = ~ mortality)


# 7. DGE DESeq2 ----
cat(paste0(">>> DGE DESeq2 \n\n"))

# Survivor vs NonSurvivor
dds$mortality <- factor(dds$mortality, levels = c("Survivor","NonSurvivor"))
dds$mortality <- relevel(dds$mortality, ref = "NonSurvivor")
dds1 <- DESeq(dds)

# NonSurvivor vs Survivor
dds$mortality <- relevel(dds$mortality, ref = "Survivor")
dds2 <- DESeq(dds)

# optional lines
dds$mortality <- relevel(dds$mortality, ref = "NonSurvivor")
res_names <- resultsNames(dds1)
res_names2 <- resultsNames(dds2)

# 8. format DGE results ----

# survivor (1) vs. non-survivor (2)
res <- results(dds1)
res_df <- res[order(res$padj),]
res_df <- as.data.frame(res_df)

res_df$gene_ensembl <- row.names(res_df)
res_df$gene         <- tstrsplit(res_df$gene_ensembl,"_")[[1]]

res_df_sig <- res_df[res_df$padj < 0.05,]

# non-survivor vs. survivor
res2 <- results(dds2)
res2_df <- res2[order(res2$padj),]
res2_df <- as.data.frame(res2_df)

res2_df$gene_ensembl <- row.names(res2_df)
res2_df$gene         <- tstrsplit(res2_df$gene_ensembl,"_")[[1]]

res2_df_sig <- res2_df[res2_df$padj < 0.05,]

# 9. BCV calculation ----
mean_counts <- rowMeans(counts(dds1, normalized = TRUE))
disp <- dispersions(dds1) 
bcv <- sqrt(disp)

df_bcv <- data.frame(
  mean_count = mean_counts,
  dispersion = disp,
  BCV = bcv
)
df_bcv$mean_count_log10 <- log10(df_bcv$mean_count)

fitted_disp <- dispersionFunction(dds1)(df_bcv$mean_count)
df_bcv$fitted_BCV <- sqrt(fitted_disp)

# Make ggplot
bcv_plot <- ggplot(df_bcv, aes(x = mean_count_log10, y = BCV)) +
  geom_point(alpha = 0.4, size = 1, color = "gray40") +
  geom_line(aes(y = fitted_BCV), color = "red", linewidth = 1.2) +
  labs(
    title = "BCV vs Mean Expression (DESeq2)",
    x = "log10(Mean normalized count)",
    y = "BCV"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )


# 10. Write CSVs ----
write.csv(res_df,file = dge_results_file1,quote = F,row.names = F)
write.csv(res_df_sig,file = dge_results_file1_sig,quote = F,row.names = F)

# non-survivor vs survivor
write.csv(res2_df,file = dge_results_file1_rev,quote = F,row.names = F)
write.csv(res2_df_sig,file = dge_results_file1_sig_rev,quote = F,row.names = F)

# BCV related things
write.csv(df_bcv,file = bcv_results_fh,row.names = F,quote = F)
pdf(file = file.path(projectdirectory,"plots_for_figures","supp_figures","BCV_plot_zi_deseq_survivor_vs_nonsurvivor.pdf"),
    height = 5,width = 6)
print(bcv_plot)
dev.off()

