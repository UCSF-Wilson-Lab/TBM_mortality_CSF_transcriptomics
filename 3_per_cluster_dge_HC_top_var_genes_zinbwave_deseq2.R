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
hc_cluster_ids_df_fh          <- file.path(projectdirectory,"tables","cluster_ID_table_hierarchical_clustering.csv")

# Output DGE 
output_dir            <- file.path(projectdirectory,"tables","DGE_results")
dir.create(output_dir)

# DeSeq2
dge_results_file1 <- file.path(output_dir,"cluster1_dge_results_zi_deseq.csv")
dge_results_file2 <- file.path(output_dir,"cluster2_dge_results_zi_deseq.csv")
dge_results_file3 <- file.path(output_dir,"cluster3_dge_results_zi_deseq.csv")


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

# Load hierarchical cluster IDs
cluster_df <- read.csv(hc_cluster_ids_df_fh,stringsAsFactors = F)

metadata   <- merge(metadata,cluster_df,by="Sample_ID")

# Make cluster DGE columns

### Cluster IDs of interest ###
cluster_column_annotations <- metadata$cluster_id_hc_top_var
###
metadata$cluster1_dge <- cluster_column_annotations
metadata$cluster2_dge <- cluster_column_annotations
metadata$cluster3_dge <- cluster_column_annotations

metadata$cluster1_dge[metadata$cluster1_dge %in% c("cluster2","cluster3")] <- "remaining"
metadata$cluster2_dge[metadata$cluster2_dge %in% c("cluster1","cluster3")] <- "remaining"
metadata$cluster3_dge[metadata$cluster3_dge %in% c("cluster2","cluster1")] <- "remaining"

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
se_obj <- SummarizedExperiment(
  assays = list(counts = df_fin), colData = data.frame(metadata), rowData = gene_meta
)

#filter out the rows that are all zeros
filter <- rowSums(assay(se_obj))
se_obj <- se_obj[filter != 0, ]


# 4. Modeling and  calculate weights ----
cat(paste0(">>> Zinbwave dimension reduction \n\n"))

MODEL_COL <- model.matrix(~cohortbatch, data=colData(se_obj) )

zinb_dim_reduc <- zinbwave(se_obj, K = 2, X=MODEL_COL, BPPARAM=BiocParallel::SerialParam(), zeroinflation = TRUE, verbose=FALSE, observationalWeights = TRUE )
weights        <- assay(zinb_dim_reduc, "weights")


# 5. Setup DGE input object and model ----

# Convert zinb dimension reduction into input objects
counts(zinb_dim_reduc)        <- as.matrix(counts(zinb_dim_reduc))

dds_c1 <- DESeqDataSet(zinb_dim_reduc, design = ~ cluster1_dge)
dds_c2 <- DESeqDataSet(zinb_dim_reduc, design = ~ cluster2_dge)
dds_c3 <- DESeqDataSet(zinb_dim_reduc, design = ~ cluster3_dge)


# Relevel and DGE analysis for each cluster
dds_c1$cluster1_dge <- relevel(dds_c1$cluster1_dge, ref = "remaining")
dds_c2$cluster2_dge <- relevel(dds_c2$cluster2_dge, ref = "remaining")
dds_c3$cluster3_dge <- relevel(dds_c3$cluster3_dge, ref = "remaining")


# 6. DGE DESeq2 ----

cat(paste0(">>> DGE DESeq2 \n\n"))

dds_c1              <- DESeq(dds_c1)
dds_c2              <- DESeq(dds_c2)
dds_c3              <- DESeq(dds_c3)


## With quantile normalization

### a. Cluster 1 ----
res1            <- results(dds_c1)
top.table.clus1 <- res1[order(res1$padj),]
top.table.clus1 <- as.data.frame(top.table.clus1)

top.table.clus1$gene_ensembl <- row.names(top.table.clus1)
top.table.clus1$gene         <- tstrsplit(top.table.clus1$gene_ensembl,"_")[[1]]
top.table.clus1$ensembl      <- tstrsplit(top.table.clus1$gene_ensembl,"_")[[2]]


### b. Cluster 2 ----
res2            <- results(dds_c2)
top.table.clus2 <- res2[order(res2$padj),]
top.table.clus2 <- as.data.frame(top.table.clus2)

top.table.clus2$gene_ensembl <- row.names(top.table.clus2)
top.table.clus2$gene         <- tstrsplit(top.table.clus2$gene_ensembl,"_")[[1]]
top.table.clus2$ensembl      <- tstrsplit(top.table.clus2$gene_ensembl,"_")[[2]]


### c. Cluster 3 ----
res3            <- results(dds_c3)
top.table.clus3 <- res3[order(res3$padj),]
top.table.clus3 <- as.data.frame(top.table.clus3)

top.table.clus3$gene_ensembl <- row.names(top.table.clus3)
top.table.clus3$gene         <- tstrsplit(top.table.clus3$gene_ensembl,"_")[[1]]
top.table.clus3$ensembl      <- tstrsplit(top.table.clus3$gene_ensembl,"_")[[2]]


# 7. Write CSVs ----
write.csv(top.table.clus1,file = dge_results_file1,quote = F,row.names = F)
write.csv(top.table.clus2,file = dge_results_file2,quote = F,row.names = F)
write.csv(top.table.clus3,file = dge_results_file3,quote = F,row.names = F)


cat(paste0("\n\n>>> DONE! <<<\n\n"))



