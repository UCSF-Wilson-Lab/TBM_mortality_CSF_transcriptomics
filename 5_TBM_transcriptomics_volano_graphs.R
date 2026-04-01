# Volcano Plot for Immune Genes in TBM Mortality
###Note: this code specifies cluster 3 as the "hypoinflammatory cluster", however this is actually annotated as cluster 1 in the manuscript (i.e. clusters 1 and 3 are switched in the manuscript)
# ============================================================================

# Load required libraries
library(EnhancedVolcano)
library(edgeR)
library(data.table)
library(dplyr)
library(gridExtra)

library(reactome.db)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Define immune-related gene keywords

# Pull Entrez IDs mapping to the Immune System pathway
path2gene <- AnnotationDbi::toTable(reactomePATHID2EXTID)
path2name <- AnnotationDbi::toTable(reactomePATHID2NAME)

# Collect pathway IDs: direct match to R-HSA-168256 plus any pathway
# whose name contains "immune" to catch nested sub-pathways
immune_pathway_ids <- path2name %>%
  filter(
    DB_ID == "R-HSA-168256" |
      grepl("immune", path_name, ignore.case = TRUE)
  ) %>%
  pull(DB_ID) %>%
  unique()


# Map Entrez IDs to HGNC symbols
sym_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = as.character(immune_entrez),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

immune_genes <- data.frame(
  entrez_id   = names(sym_map),
  gene_symbol = unname(sym_map),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct() %>%
  arrange(gene_symbol)

cat("reactome.db version :", as.character(packageVersion("reactome.db")), "\n")
cat("Immune genes found  :", nrow(immune_genes), "\n")


# Create combined immune gene pattern
immune_pattern <- paste(immune_genes$gene_symbol, collapse = "|")

immune_pattern
# ============================================================================
# Set up directories
# ============================================================================
projectdirectory      <- "/data/" #your path here
resouce_dir           <- file.path(projectdirectory, "resources")
input_dir             <- file.path(projectdirectory, "input")
plot_dir              <- file.path(projectdirectory, "plots_for_figures", "main_figures")
plot_dir_supp         <- file.path(projectdirectory, "plots_for_figures", "supp_figures")
table_dir             <- file.path(projectdirectory, "tables")
object_dir            <- file.path(projectdirectory, "objects")

# Input files
processed_gene_counts_file <- file.path(input_dir, "input_gene_count_matrix.csv")
metadata_tbm_file          <- file.path(input_dir, "metadata_TBM_formatted.csv")
biomart_metadata_fh         <- file.path(resouce_dir,"biomart_protein_coding_genes.tsv")

# DGE results
dge_results_file <- file.path(table_dir, "DGE_results", "DGE_results_ZI_batch_opt_model", "dge_results_zi_deseq.csv")

# Output directories
output <- '/data/sreddy/TBM_transcriptomic_figures'
output_dir <- file.path(output, "volcano", "immune_genes")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load and Process Data
# ============================================================================
df <- read.csv(processed_gene_counts_file, row.names = 1, check.names = F)
gene_metadata <- read.delim(biomart_metadata_fh, stringsAsFactors = F)
metadata <- read.csv(metadata_tbm_file, stringsAsFactors = FALSE)
row.names(metadata) <- metadata$Sample_ID
metadata <- metadata[metadata$Group_rx %in% "1",]
metadata <- metadata[metadata$Sample_ID %in% names(df),]
metadata$Mortality <- as.factor(metadata$short_term_Mortality)
metadata$Gender <- as.factor(metadata$gender)
metadata$Batch <- as.numeric(metadata$BatchAnalysis)
metadata$cohortbatch <- paste("TBM", metadata$Batch, sep = "_") %>% as.factor()
metadata$Status <- metadata$death_shortterm
metadata$Age <- metadata$age
metadata$TBM_status <- metadata$post_MNGS_TBM_status
metadata$TBM_status[metadata$TBM_status %in% "Definite TBM"] <- "DefiniteTBM"
metadata$TBM_status[metadata$TBM_status %in% "Probable TBM"] <- "ProbableTBM"
metadata$mortality <- ifelse(as.numeric(metadata$Mortality) > 1, "Non-Survivor", "Survivor")
metadata_filt <- metadata[, c("Age", "Gender", "mortality", "Status")]

# Load DGE results
res_df_sig_zi_full <- read.csv(dge_results_file, stringsAsFactors = F)
res_df_sig_zi_full$gene_name <- tstrsplit(res_df_sig_zi_full$gene_ensembl, "_")[[1]]
res_df_sig_zi <- res_df_sig_zi_full[res_df_sig_zi_full$padj < 0.05,]
res_df_sig_zi_filt <- res_df_sig_zi[res_df_sig_zi$gene_name %in% gene_metadata$hgnc_symbol,]

# ============================================================================
# Identify Immune Genes
# ============================================================================

# Identify immune genes in your dataset
res_df_sig_zi_full$is_immune <- res_df_sig_zi_full$gene_name %in% immune_genes$gene_symbol

# Subset to immune genes only
res_df_immune <- res_df_sig_zi_full[res_df_sig_zi_full$is_immune, ]

# Print summary statistics
cat("\n=== Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(res_df_sig_zi_full), "\n")
cat("Total immune genes identified:", nrow(res_df_immune), "\n")
cat("Significant immune genes (padj < 0.05):", 
    sum(res_df_immune$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated immune genes (padj < 0.05, log2FC > 1.5):", 
    sum(res_df_immune$padj < 0.05 & res_df_immune$log2FoldChange > 1.5, na.rm = TRUE), "\n")
cat("Downregulated immune genes (padj < 0.05, log2FC < -1.5):", 
    sum(res_df_immune$padj < 0.05 & res_df_immune$log2FoldChange < -1.5, na.rm = TRUE), "\n\n")


# Save immune gene list
write.csv(res_df_immune, 
          file = file.path(output_dir, "immune_genes_dge_results.csv"),
          row.names = FALSE, quote = FALSE)

# Subset to significant immune genes only (padj < 0.05, |log2FC| > 1.5)
res_df_immune_sig <- res_df_immune[
  !is.na(res_df_immune$padj) & 
    res_df_immune$padj < 0.05 & 
    abs(res_df_immune$log2FoldChange) > 1.5, 
]

up_genes <- res_df_immune_sig$gene_name[res_df_immune_sig$log2FoldChange > 1.5]
down_genes <- res_df_immune_sig$gene_name[res_df_immune_sig$log2FoldChange < -1.5]
cat("\n=== Significant Immune Genes ===\n")
cat("\nUpregulated (padj < 0.05, log2FC > 1.5):\n")
cat(paste(sort(up_genes), collapse = ", "), "\n")

cat("\nDownregulated (padj < 0.05, log2FC < -1.5):\n")
cat(paste(sort(down_genes), collapse = ", "), "\n")

# Save significant immune gene list
write.csv(res_df_immune_sig, 
          file = file.path(output_dir, "immune_genes_dge_results_significant.csv"),
          row.names = FALSE, quote = FALSE)

# ============================================================================
# Create Volcano Plots
# ============================================================================
# Define genes to label
genes_to_label <- c('FYN', 'LCK', 'LAT', 'IL7', 'CD27', 'IL12RB1', 'STAT4', 
                    'CD81', 'TNFRSF13C', 'PLCG2', 'KLRD1', 'ULBP1', 'ULBP3', 'CD8B',
                    'IL23R', 'CCR6', 'IL1A', 'TNF', 'UBE2N', 'TRAF2', 'STING1', 
                    'IRF6', 'TYK2', 'IFIT3', 'CLEC7A', 'CLEC6A', 
                    'CLEC4C', 'COLEC12', 'AIM2', 'PYCARD', 'IL1RN', 'CXCL1', 'C4A', 'C4B',
                    'UBE2L3', 'UBE2D4', 'UBE3B', 'UBE3C', 'TRIM26', 'TRIM45', 'SIAH2',
                    'PSMC2', 'PSMC3', 'PSMD2', 'CD40'
                    )

# Create capped version of ALL genes
res_df_all_capped <- res_df_sig_zi_full

# Cap p-values at -log10(p) = 5
res_df_all_capped$padj_capped <- ifelse(
  res_df_all_capped$padj < 10^-5, 
  10^-5,  
  res_df_all_capped$padj
)

# Cap log2FC at ±10
res_df_all_capped$log2FoldChange_capped <- ifelse(
  res_df_all_capped$log2FoldChange > 10,
  10,
  ifelse(
    res_df_all_capped$log2FoldChange < -10,
    -10,
    res_df_all_capped$log2FoldChange
  )
)

# Main volcano plot
plot_volcano_all_genes <- EnhancedVolcano(
  toptable = res_df_all_capped,
  lab = res_df_all_capped$gene_name,
  x = 'log2FoldChange_capped',
  y = 'padj_capped',
  ylim = c(0, 5.0),  # Limit y-axis display to 0-7.5
  title = 'Differentially Expressed Genes: TBM Survivors vs Nonsurvivors',
  subtitle = expression("Adjusted p-values capped at -log"[10]*"(p)=5.0, log"[2]*"FC capped at ±10 for clarity"),
  selectLab = genes_to_label,
  pCutoff = 0.05,
  FCcutoff = 1.5,
  pointSize = 2.0,
  labSize = 5.5,
  colAlpha = 0.5,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  boxedLabels = TRUE,
  max.overlaps = Inf,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  labCol = 'black',
  labFace = 'bold'
)

plot_volcano_all_genes <- plot_volcano_all_genes +
  ylab(expression("-Log"[10]*" (adjusted p-value)"))

plot_volcano_all_genes <- plot_volcano_all_genes +
  scale_color_manual(
    values = c("grey30", "forestgreen", "royalblue", "red2"),
    labels = c("Not significant", "Log2FC > |1.5|", "Adj p-value < 0.05", "Adj p-value < 0.05 & Log2FC > |1.5|")
  )

plot_volcano_all_genes <- plot_volcano_all_genes +
  theme(
    legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )


plot_volcano_all_genes <- plot_volcano_all_genes + 
  theme(
    plot.title = element_text(size = 25, face = "bold"),          # Title
    plot.subtitle = element_text(size = 18, face = "italic"),     # Subtitle
    axis.title = element_text(size = 1, face = "bold"),          # Both axis labels
    axis.title.x = element_text(size = 16, face = "bold"),        # X-axis label only
    axis.title.y = element_text(size = 16, face = "bold"),        # Y-axis label only
    axis.text = element_text(size = 16)                           # Axis tick labels
  )

pdf(file = file.path(output_dir, "volcano_all_genes_immune_labeled.pdf"), 
    height = 10, width = 12)
print(plot_volcano_all_genes)
dev.off()

#===================================
# Volcano plot, cluster 3 vs clusters 1&2
# Cluster 3
dge_hclust_cluster_results_fh3 <- file.path(table_dir,"DGE_results","patient_Hclust_top_var_clusters_DGE_results_ZI_batch_opt_model","cluster3_dge_results_zi_deseq.csv")
top.table.hclus3.full <- read.csv(dge_hclust_cluster_results_fh3, stringsAsFactors = F)
top.table.hclus3      <- top.table.hclus3.full[top.table.hclus3.full$padj < 0.05,]
top.table.hclus3$absFC <- abs(top.table.hclus3$log2FoldChange)
top.table.hclus3       <- top.table.hclus3[top.table.hclus3$absFC >= 2,]
top.table.hclus3 <- top.table.hclus3[order(top.table.hclus3$log2FoldChange, decreasing = T),]
top.table.hclus3$cluster_id <- "cluster3"

# Identify immune genes in your dataset
top.table.hclus3.full$is_immune <- top.table.hclus3.full$gene %in% immune_genes$gene_symbol


# Subset to immune genes only
cluster3_immune <- top.table.hclus3.full[top.table.hclus3.full$is_immune, ]

# Print summary statistics
cat("\n=== Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(top.table.hclus3.full), "\n")
cat("Total immune genes identified:", nrow(cluster3_immune), "\n")
cat("Significant immune genes (padj < 0.05):", 
    sum(cluster3_immune$padj < 0.05, na.rm = TRUE), "\n")


# Save immune gene list
write.csv(cluster3_immune, 
          file = file.path(output_dir, "cluster_immune_genes_dge_results.csv"),
          row.names = FALSE, quote = FALSE)

# Subset to significant immune genes only (padj < 0.05, |log2FC| > 1.5)
cluster3_immune_sig <- cluster3_immune[
  !is.na(cluster3_immune$padj) & 
    cluster3_immune$padj < 0.05 & 
    abs(cluster3_immune$log2FoldChange) > 1.5, 
]

up_genes <- cluster3_immune_sig$gene[cluster3_immune_sig$log2FoldChange > 1.5]
down_genes <- cluster3_immune_sig$gene[cluster3_immune_sig$log2FoldChange < -1.5]
cat("\n=== Significant Immune Genes ===\n")
cat("\nUpregulated (padj < 0.05, log2FC > 1.5):\n")
cat(paste(sort(up_genes), collapse = ", "), "\n")
cat("\nDownregulated (padj < 0.05, log2FC < -1.5):\n")
cat(paste(sort(down_genes), collapse = ", "), "\n")

# Save significant immune gene list
write.csv(res_df_immune_sig, 
          file = file.path(output_dir, "cluster_immune_genes_dge_results_significant.csv"),
          row.names = FALSE, quote = FALSE)

# Create capped version
top.table.hclus3.capped <- top.table.hclus3.full

# Cap p-values at -log10(p) = 10, which corresponds to p = 10^-10
top.table.hclus3.capped$padj_capped <- ifelse(
  top.table.hclus3.capped$padj < 10^-10, 
  10^-10,
  top.table.hclus3.capped$padj
)

# Cap log2FC at ±10
top.table.hclus3.capped$log2FoldChange_capped <- ifelse(
  top.table.hclus3.capped$log2FoldChange > 10,
  20,
  ifelse(
    top.table.hclus3.capped$log2FoldChange < -10,
    -10,
    top.table.hclus3.capped$log2FoldChange
  )
)

cluster3_genes_label <- c('HLA-DRA', 'HLA-DRB1', 'TAP1', 'B2M', 'LAT', 'CD3G', 
                          'ITGAL', 'IL12RB1', 'STAT1', 'IFNGR2', 'CD14', 'CSF1R',
                          'TLR4', 'C1QA', 'C1QC', 'C7', 'CXCL2', 'FCGR3B', 
                          'CD14', 'CD40', 'CTLA4', 'IL27RA', 'CCR1', 'CCR6', 
                          'FCGR1A', 'ITGA4', 'STAT5B', 'IL17RA', 'IL1B', 'REL', 
                          "B2M", 'TAP1', 'C1QA', 'C1QC', 'C7', 'CCL3', 
                          'TGFBR1', 'CD8A', 'IFNAR1', 'HLA-E', 'CXCL2', 'FCGR3B', 'CD16b')

plot_volcano_zi_hc3 <- EnhancedVolcano(
  toptable = top.table.hclus3.capped,
  lab = top.table.hclus3.capped$gene,
  x = 'log2FoldChange_capped',
  y = 'padj_capped',
  ylim = c(0, 10),  # Limit y-axis to 0-10
  title = 'Hypoinflammatory Cluster',
  subtitle = 'Differentially expressed genes',
  selectLab = cluster3_genes_label,
  pCutoff = 0.05,
  FCcutoff = 1.5,
  pointSize = 2.0,
  labSize = 4.0,
  colAlpha = 0.5,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  boxedLabels = TRUE,
  max.overlaps = Inf,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  labCol = 'black',
  labFace = 'bold',
)


plot_volcano_zi_hc3  <- plot_volcano_zi_hc3 +
  ylab(expression("-Log"[10]*" (adjusted p-value)"))


plot_volcano_zi_hc3 <- plot_volcano_zi_hc3 + 
  theme(
    plot.title = element_text(size = 25, face = "bold"),          # Title
    plot.subtitle = element_text(size = 16, face = "italic"),     # Subtitle
    axis.title = element_text(size = 1, face = "bold"),          # Both axis labels
    axis.title.x = element_text(size = 16, face = "bold"),        # X-axis label only
    axis.title.y = element_text(size = 16, face = "bold"),        # Y-axis label only
    axis.text = element_text(size = 16)                           # Axis tick labels
  )

pdf(file = file.path(output_dir, "fig3b_volcano_tbm_hclust_cluster3.pdf"), 
    height = 8, width = 8)
print(plot_volcano_zi_hc3)
dev.off()


#===============volcano plot for cluster 3 mortality. 
top.table.hclust3.survival <- read.csv('/data/sreddy/cluster3/tables/DGE_results/survivors_cluster3_DGE_results_ZI_batch_opt_model/survivors_cluster3_dge_results_zi_deseq.csv', stringsAsFactors = F)


# Identify immune genes in your dataset
top.table.hclust3.survival$is_immune <- top.table.hclust3.survival$gene %in% immune_genes$gene_symbol

# Subset to immune genes only
cluster3surv_immune <- top.table.hclust3.survival[top.table.hclust3.survival$is_immune, ]
# Print summary statistics
cat("\n=== Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(top.table.hclust3.survival), "\n")
cat("Total immune genes identified:", nrow(cluster3surv_immune), "\n")
cat("Significant immune genes (padj < 0.05):", 
    sum(cluster3surv_immune$padj < 0.05, na.rm = TRUE), "\n")

sig_immune_only_surv <- cluster3surv_immune$gene[
  cluster3_immune$padj < 0.05 & !is.na(cluster3_immune$padj)
]

# Save immune gene list
write.csv(cluster3surv_immune, 
          file = file.path(output_dir, "cluster_surv_immune_genes_dge_results.csv"),
          row.names = FALSE, quote = FALSE)

# Subset to significant immune genes only (padj < 0.05, |log2FC| > 1.5)
cluster3surv_immune_sig <- cluster3surv_immune[
  !is.na(cluster3surv_immune$padj) & 
    cluster3surv_immune$padj < 0.05 & 
    abs(cluster3surv_immune$log2FoldChange) > 1.5, 
]

up_genes <- cluster3surv_immune_sig$gene[cluster3surv_immune_sig$log2FoldChange > 1.5]
down_genes <- cluster3surv_immune_sig$gene[cluster3surv_immune_sig$log2FoldChange < -1.5]
cat("\n=== Significant Immune Genes ===\n")
cat("\nUpregulated (padj < 0.05, log2FC > 1.5):\n")
cat(paste(sort(up_genes), collapse = ", "), "\n")
cat("\nDownregulated (padj < 0.05, log2FC < -1.5):\n")
cat(paste(sort(down_genes), collapse = ", "), "\n")

# Save significant immune gene list
write.csv(res_df_immune_sig, 
          file = file.path(output_dir, "cluster_surv_immune_genes_dge_results_significant.csv"),
          row.names = FALSE, quote = FALSE)


# Create capped version
top.table.hclust3.survival.capped <- top.table.hclust3.survival

# Cap p-values at -log10(p) = 5,
top.table.hclust3.survival.capped$padj_capped <- ifelse(
  top.table.hclust3.survival.capped$padj < 10^-5, 
  10^-5,
  top.table.hclust3.survival.capped$padj
)

# Cap log2FC at ±10
top.table.hclust3.survival.capped$log2FoldChange_capped <- ifelse(
  top.table.hclust3.survival.capped$log2FoldChange > 10,
  10,
  ifelse(
    top.table.hclust3.survival.capped$log2FoldChange < -10,
    -10,
    top.table.hclust3.survival.capped$log2FoldChange
  )
)

genes_lab <- c('TNFRSF1A', 'TNFRSF1B', 'IFNGR2', 'ITGAL',
               'FCER1G', 'CCL5', 'IL1RAP', 'TNFAIP8', 'HLA-F', 'KLRD1', 'PAX5',
               'CXCL16', 'HLA-DPB1')

plot_volcano_hc3_surv <- EnhancedVolcano(
  toptable = top.table.hclust3.survival.capped,
  lab = top.table.hclust3.survival.capped$gene,
  x = 'log2FoldChange_capped',
  y = 'padj_capped',
  ylim = c(0, 5),  # Limit y-axis to 0-10
  title = 'HC: TBM Survivors vs Nonsurvivors',
  subtitle = 'Differentially expressed genes in comparison to cluster nonsurvivors ',
  selectLab = genes_lab,
  pCutoff = 0.05,
  FCcutoff = 1.5,
  pointSize = 2.0,
  labSize = 4.0,
  colAlpha = 0.5,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  boxedLabels = TRUE,
  max.overlaps = Inf,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  labCol = 'black',
  labFace = 'bold',
)

plot_volcano_hc3_surv  <- plot_volcano_hc3_surv +
  ylab(expression("-Log"[10]*" (adjusted p-value)"))

plot_volcano_hc3_surv <- plot_volcano_hc3_surv + 
  theme(
    plot.title = element_text(size = 25, face = "bold"),          # Title
    plot.subtitle = element_text(size = 16, face = "italic"),     # Subtitle
    axis.title = element_text(size = 1, face = "bold"),          # Both axis labels
    axis.title.x = element_text(size = 16, face = "bold"),        # X-axis label only
    axis.title.y = element_text(size = 16, face = "bold"),        # Y-axis label only
    axis.text = element_text(size = 16)                           # Axis tick labels
  )

pdf(file = file.path(output_dir, "cluster3_survival.pdf"), 
    height = 8, width = 8)
print(plot_volcano_hc3_surv)
dev.off()

#===============================================================================
# Volcano plot: Survivors vs Non-Survivors in Definite TBM cases only
#===============================================================================

# Load DGE results for definite TBM only
dge_definite_only <- file.path(table_dir, "DGE_results", "DefTBM_DGE_results_ZI_batch_opt_model", "dge_results_zi_deseq.csv")
top.table.definite.full <- read.csv(dge_definite_only, stringsAsFactors = F)

# Identify immune genes
top.table.definite.full$is_immune <- top.table.definite.full$gene %in% immune_genes$gene_symbol

# Subset to immune genes only
definite_immune <- top.table.definite.full[top.table.definite.full$is_immune, ]

# Print summary statistics
cat("\n=== Definite TBM Only - Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(top.table.definite.full), "\n")
cat("Total immune genes identified:", nrow(definite_immune), "\n")
cat("Significant immune genes (padj < 0.05):",
    sum(definite_immune$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated immune genes (padj < 0.05, log2FC > 1.5):",
    sum(definite_immune$padj < 0.05 & definite_immune$log2FoldChange > 1.5, na.rm = TRUE), "\n")
cat("Downregulated immune genes (padj < 0.05, log2FC < -1.5):",
    sum(definite_immune$padj < 0.05 & definite_immune$log2FoldChange < -1.5, na.rm = TRUE), "\n\n")

# Save immune gene list
write.csv(definite_immune,
          file = file.path(output_dir, "definite_tbm_immune_genes_dge_results.csv"),
          row.names = FALSE, quote = FALSE)

# Subset to significant immune genes only (padj < 0.05, |log2FC| > 1.5)
definite_immune_sig <- definite_immune[
  !is.na(definite_immune$padj) & 
    definite_immune$padj < 0.05 & 
    abs(definite_immune$log2FoldChange) > 1.5, 
]

up_genes <- definite_immune_sig$gene[definite_immune_sig$log2FoldChange > 1.5]
down_genes <- definite_immune_sig$gene[definite_immune_sig$log2FoldChange < -1.5]
cat("\n=== Significant Immune Genes ===\n")
cat("\nUpregulated (padj < 0.05, log2FC > 1.5):\n")
cat(paste(sort(up_genes), collapse = ", "), "\n")

cat("\nDownregulated (padj < 0.05, log2FC < -1.5):\n")
cat(paste(sort(down_genes), collapse = ", "), "\n")

# Save significant immune gene list
write.csv(res_df_immune_sig, 
          file = file.path(output_dir, "immune_genes_dge_results_significant.csv"),
          row.names = FALSE, quote = FALSE)

# Create capped version
top.table.definite.capped <- top.table.definite.full

# Cap p-values at -log10(p) = 5
top.table.definite.capped$padj_capped <- ifelse(
  top.table.definite.capped$padj < 10^-5,
  10^-5,
  top.table.definite.capped$padj
)

# Cap log2FC at ±10
top.table.definite.capped$log2FoldChange_capped <- ifelse(
  top.table.definite.capped$log2FoldChange > 10,
  10,
  ifelse(
    top.table.definite.capped$log2FoldChange < -10,
    -10,
    top.table.definite.capped$log2FoldChange
  )
)

genes_to_label_definite <- c('FYN', 'LCK', 'IL7', 'IL12RB1', 'STAT4', 
                    'CD8B','IL23R', 'CCR6', 'TRAF2', 'STING1',
                    'IL1RN', 'CXCL1', 'C4A', 'C4B',
                    'UBE2D4', 'UBE3C', 'TRIM26', 'SIAH2',
                    'CD40', 'C1QB', 'C1QC', 'C2', 'CFB', 
                    'KLRK1', 'IRF7', 'IFI27', 'STING1'
)


# Generate volcano plot
plot_volcano_definite <- EnhancedVolcano(
  toptable = top.table.definite.capped,
  lab = top.table.definite.capped$gene,
  x = 'log2FoldChange_capped',
  y = 'padj_capped',
  ylim = c(0, 5.0),
  title = 'Definite TBM: Survivors vs Non-Survivors',
  subtitle = 'Adjusted p-values capped at -log10(p)=5.0, log2FC capped at ±10 for clarity',
  selectLab = genes_to_label_definite,
  pCutoff = 0.05,
  FCcutoff = 1.5,
  pointSize = 2.0,
  labSize = 5.5,
  colAlpha = 0.5,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  boxedLabels = TRUE,
  max.overlaps = Inf,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  labCol = 'black',
  labFace = 'bold',
)

plot_volcano_definite <- plot_volcano_definite +
  ylab(expression("-Log"[10]*" (adjusted p-value)"))

plot_volcano_definite <- plot_volcano_definite +
  theme(
    plot.title = element_text(size = 25, face = "bold"),
    plot.subtitle = element_text(size = 18, face = "italic"),
    axis.title = element_text(size = 1, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16)
  )

pdf(file = file.path(output_dir, "volcano_definite_tbm_survivors_immune.pdf"),
    height = 10, width = 12)
print(plot_volcano_definite)
dev.off()



#===============================================================================
# Volcano plot: Definite vs Probable TBM
#===============================================================================
# Load DGE results for TBM status comparison (ProbableTBM vs DefiniteTBM)
top.table.tbm_status.full <- read.csv('~/TBM_status/tables/DGE_results/TBM_status_DGE_results_ZI_batch_opt_model/TBM_status_dge_results_zi_deseq.csv', stringsAsFactors = F)

# Identify immune genes
top.table.tbm_status.full$is_immune <- top.table.tbm_status.full$gene %in% immune_genes$gene_symbol

# Subset to immune genes only
tbm_status_immune <- top.table.tbm_status.full[top.table.tbm_status.full$is_immune, ]

# Print summary statistics
cat("\n=== TBM Status (ProbableTBM vs DefiniteTBM) - Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(top.table.tbm_status.full), "\n")
cat("Total immune genes identified:", nrow(tbm_status_immune), "\n")
cat("Significant immune genes (padj < 0.05):",
    sum(tbm_status_immune$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated immune genes (padj < 0.05, log2FC > 1.5):",
    sum(tbm_status_immune$padj < 0.05 & tbm_status_immune$log2FoldChange > 1.5, na.rm = TRUE), "\n")
cat("Downregulated immune genes (padj < 0.05, log2FC < -1.5):",
    sum(tbm_status_immune$padj < 0.05 & tbm_status_immune$log2FoldChange < -1.5, na.rm = TRUE), "\n\n")

# Save immune gene list
write.csv(tbm_status_immune,
          file = file.path(output_dir, "tbm_status_immune_genes_dge_results.csv"),
          row.names = FALSE, quote = FALSE)

# Subset to significant immune genes only (padj < 0.05, |log2FC| > 1.5)
tbm_status_immune_sig <- tbm_status_immune[
  !is.na(tbm_status_immune$padj) &
    tbm_status_immune$padj < 0.05 &
    abs(tbm_status_immune$log2FoldChange) > 1.5,
]

up_genes   <- tbm_status_immune_sig$gene[tbm_status_immune_sig$log2FoldChange > 1.5]
down_genes <- tbm_status_immune_sig$gene[tbm_status_immune_sig$log2FoldChange < -1.5]

cat("\n=== Significant Immune Genes ===\n")
cat("\nUpregulated in ProbableTBM (padj < 0.05, log2FC > 1.5):\n")
cat(paste(sort(up_genes), collapse = ", "), "\n")
cat("\nDownregulated in ProbableTBM (padj < 0.05, log2FC < -1.5):\n")
cat(paste(sort(down_genes), collapse = ", "), "\n")

# Save significant immune gene list
write.csv(tbm_status_immune_sig,
          file = file.path(output_dir, "tbm_status_immune_genes_dge_results_significant.csv"),
          row.names = FALSE, quote = FALSE)

# Create capped version
top.table.tbm_status.capped <- top.table.tbm_status.full

# Cap p-values at -log10(p) = 5
top.table.tbm_status.capped$padj_capped <- ifelse(
  top.table.tbm_status.capped$padj < 10^-5,
  10^-5,
  top.table.tbm_status.capped$padj
)

# Cap log2FC at ±10
top.table.tbm_status.capped$log2FoldChange_capped <- ifelse(
  top.table.tbm_status.capped$log2FoldChange > 10,
  10,
  ifelse(
    top.table.tbm_status.capped$log2FoldChange < -10,
    -10,
    top.table.tbm_status.capped$log2FoldChange
  )
)

gene_list_def_prob <- c('IL1A','IL1RAP', 'CLEC4A', 'CLEC4D',
                        'CLEC4E', 'CLEC5A', 'S100A8', 'S100A9',
                        'S100A12', 'CXCL8', 'CCL2', 'CCL3'
                        )

# Generate volcano plot
plot_volcano_tbm_status <- EnhancedVolcano(
  toptable  = top.table.tbm_status.capped,
  lab       = top.table.tbm_status.capped$gene,
  x         = 'log2FoldChange_capped',
  y         = 'padj_capped',
  ylim      = c(0, 5.0),
  title     = 'TBM Status: Probable TBM vs Definite TBM (Select genes)',
  subtitle  = 'Adjusted p-values capped at -log10(p)=5.0, log2FC capped at ±10 for clarity',
  selectLab = gene_list_def_prob,
  pCutoff   = 0.05,
  FCcutoff  = 1.5,
  pointSize = 2.0,
  labSize   = 5.5,
  colAlpha  = 0.5,
  drawConnectors  = TRUE,
  widthConnectors = 0.5,
  colConnectors   = 'black',
  boxedLabels     = TRUE,
  max.overlaps    = Inf,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  labCol  = 'black',
  labFace = 'bold',
)

plot_volcano_tbm_status <- plot_volcano_tbm_status +
  ylab(expression("-Log"[10]*" (adjusted p-value)"))

plot_volcano_tbm_status <- plot_volcano_tbm_status +
  theme(
    plot.title    = element_text(size = 25, face = "bold"),
    plot.subtitle = element_text(size = 18, face = "italic"),
    axis.title    = element_text(size = 1,  face = "bold"),
    axis.title.x  = element_text(size = 16, face = "bold"),
    axis.title.y  = element_text(size = 16, face = "bold"),
    axis.text     = element_text(size = 16)
  )

pdf(file   = file.path(output_dir, "volcano_tbm_status_immune.pdf"),
    height = 10, width = 12)
print(plot_volcano_tbm_status)
dev.off()



