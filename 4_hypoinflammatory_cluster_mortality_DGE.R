###Note: this code specifies cluster 3 as the "hypoinflammatory cluster", however this is actually annotated as cluster 1 in the manuscript (i.e. clusters 1 and 3 are switched in the manuscript)

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
projectdirectory  <- "/data/rdandekar/rprojects/TBM_RNASeq_Analysis"
input_dir         <- file.path(projectdirectory,"input")
resource_dir      <- file.path(projectdirectory,"resources")
processed_gene_counts_file    <- file.path(input_dir,"input_gene_count_matrix.csv")
metadata_tbm_file             <- file.path(input_dir, "metadata_TBM_formatted.csv")
biomart_metadata_fh           <- file.path(resource_dir,"biomart_protein_coding_genes.tsv")
cluster_coord_df_fh           <- file.path(projectdirectory,"tables","coords_sample_UMAP_all_genes.csv")
hc_cluster_ids_df_fh          <- file.path(projectdirectory,"tables","cluster_ID_table_hierarchical_clustering_Oct2025.csv")

# Output DGE 
output <- '/data/sreddy/cluster3'
output_dir            <- file.path(output,"tables","DGE_results","survivors_cluster3_DGE_results_ZI_batch_opt_model")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# DeSeq2
dge_results_file <- file.path(output_dir,"survivors_cluster3_dge_results_zi_deseq.csv")


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

# Add in cluster IDs
cluster_df <- read.csv(cluster_coord_df_fh,stringsAsFactors = F)
cluster_df <- cluster_df[,c("Sample_ID","cluster_id","V1","V2","V3")]
# Load hierarchical cluster IDs
hc_cluster_df <- read.csv(hc_cluster_ids_df_fh,stringsAsFactors = F)
# merge all cluster IDs into one dataframe
cluster_df <- merge(cluster_df,hc_cluster_df,by="Sample_ID")

metadata <- merge(metadata,cluster_df,by="Sample_ID")

# Statistical comparison of nonzero protein coding genes across clusters
cat(paste0(">>> Statistical Analysis: Nonzero Protein Coding Genes Across Clusters\n\n"))

# Calculate summary statistics by cluster (before filtering to cluster 3)
summary_by_cluster <- metadata %>%
  group_by(cluster_id_hc_top_var) %>%
  summarise(
    n = n(),
    mean_nonzero = mean(numNONZero),
    median_nonzero = median(numNONZero),
    sd_nonzero = sd(numNONZero),
    q25_nonzero = quantile(numNONZero, 0.25),
    q75_nonzero = quantile(numNONZero, 0.75),
    IQR_nonzero = IQR(numNONZero)
  )

print(summary_by_cluster, width = Inf)

# ANOVA
anova_result <- aov(numNONZero ~ cluster_id_hc_top_var, data = metadata)
anova_summary <- summary(anova_result)
cat(paste0("ANOVA p-value: ", format(anova_summary[[1]][["Pr(>F)"]][1], scientific = TRUE), "\n\n"))


# FILTER FOR CLUSTER 3 ONLY
cat(paste0(">>> Filtering for cluster 3 samples only \n\n"))
metadata <- metadata[metadata$cluster_id_hc_top_var == "cluster3", ]
cat(paste0("Number of cluster 3 samples: ", nrow(metadata), "\n"))
cat(paste0("Survivors: ", sum(metadata$mortality == "Survivor"), "\n"))
cat(paste0("Non-survivors: ", sum(metadata$mortality == "NonSurvivor"), "\n\n"))

# Convert mortality to factor for DESeq2
metadata$mortality <- as.factor(metadata$mortality)

# 2. Filter out genes ----
# * Minimum QC is to require non-zero counts in at least 20% of samples per gene
# * If this is too strict, then we can set the CUTOFF to 1
CUTOFF = 0.2*nrow(metadata)  # Changed from ncol(df) to reflect cluster 3 sample count

# Filter gene counts to only include cluster 3 samples
df_cluster3 <- df[, metadata$Sample_ID]

# Filter for genes detected in >=20% of cluster 3 samples
keep <- rowSums(df_cluster3>0) >= CUTOFF
keep <- names(keep[keep == TRUE])
total_genes_kept <- length(keep) # number of genes retained
cat(paste0("Genes kept after filtering: ", total_genes_kept, "\n\n"))

# Filter DGE object
df_filt <- df_cluster3[keep, ]
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

#filter out the rows that are all zeros
filter <- rowSums(assay(se_obj))
se_obj <- se_obj[filter != 0, ]


# 4. Add Quartile categories ----

quartiles <- quantile(se_obj$numNONZero, probs = c(0, 0.25, 0.5, 0.75, 1))

categorize_quartile <- function(value, quartiles) {
  if (value <= quartiles[2]) {
    return("Q1")
  } else if (value <= quartiles[3]) {
    return("Q2")
  } else if (value <= quartiles[4]) {
    return("Q3")
  } else {
    return("Q4")
  }
}

# Apply the function to create a new variable
se_obj$Quartilenonzero <- sapply(se_obj$numNONZero, categorize_quartile, quartiles = quartiles)


# 5. Modeling and calculate weights ----
cat(paste0(">>> Zinbwave dimension reduction \n\n"))

# Optimal Model for cluster 3 mortality analysis
MODEL_COL <- model.matrix(~cohortbatch, data=colData(se_obj))
MODEL_ROW <- model.matrix(~gccontent + log(length), data = rowData(se_obj))

zinb_dim_reduc <- zinbwave(se_obj, K = 2, X=MODEL_COL, BPPARAM=BiocParallel::SerialParam(), 
                           zeroinflation = TRUE, verbose=FALSE, observationalWeights = TRUE)
weights        <- assay(zinb_dim_reduc, "weights")


# 6. Setup DGE input object and model ----

# Convert zinb dimension reduction into input objects
counts(zinb_dim_reduc) <- as.matrix(counts(zinb_dim_reduc))

# Create DESeq2 dataset with mortality as the design variable
dds <- DESeqDataSet(zinb_dim_reduc, design = ~ mortality)

# Relevel to make "NonSurvivor" the reference group
dds$mortality <- relevel(dds$mortality, ref = "NonSurvivor")


# 7. DGE DESeq2 ----

cat(paste0(">>> DGE DESeq2 - Mortality comparison within cluster 3\n\n"))

dds <- DESeq(dds)


# 8. Extract results ----

res <- results(dds, contrast = c("mortality", "Survivor", "NonSurvivor"))
top.table <- res[order(res$padj), ]
top.table <- as.data.frame(top.table)

top.table$gene_ensembl <- row.names(top.table)
top.table$gene         <- tstrsplit(top.table$gene_ensembl,"_")[[1]]
top.table$ensembl      <- tstrsplit(top.table$gene_ensembl,"_")[[2]]

# Add summary statistics
cat(paste0("Significant genes (padj < 0.05): ", sum(top.table$padj < 0.05, na.rm = TRUE), "\n"))
cat(paste0("Significant genes (padj < 0.01): ", sum(top.table$padj < 0.01, na.rm = TRUE), "\n\n"))




# 9. Write CSVs ----
write.csv(top.table, file = dge_results_file, quote = F, row.names = F)

# Weight CSV
write.csv(weights, file = file.path(output_dir,"matrix_weights_zinb_ALL_GENES.csv"), row.names = T, quote = F)

cat(paste0("\n\n>>> DONE! <<<\n\n"))
