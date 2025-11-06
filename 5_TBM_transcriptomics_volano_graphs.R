# Volcano Plot for Immune Genes in TBM Mortality
###Note: this code specifies cluster 3 as the "hypoinflammatory cluster", however this is actually annotated as cluster 1 in the manuscript (i.e. clusters 1 and 3 are switched in the manuscript)
# ============================================================================

# Load required libraries
library(EnhancedVolcano)
library(edgeR)
library(data.table)
library(dplyr)
library(gridExtra)

# Define immune-related gene keywords
immune_keywords <- list(
  # Cytokines and interleukins
  cytokines = c("^IL[0-9]+[A-Z]*$", "^TNF", "^TNFA", "^TNFB",
                "^IFN", "^IFNA", "^IFNB", "^IFNG", 
                "^CSF[0-9]", "^GM-CSF", "^M-CSF", "^LIF$", "^OSM$"),
  
  # Cytokine receptors
  cytokine_receptors = c("^IL[0-9]+R[A-Z]*$", "^IFNGR", "^IFNAR", "^IFNLR",
                         "^TNFRSF", "^TNFAIP", "^CSF[0-9]R", "^TGFBR"),
  
  # Chemokines
  chemokines = c("^CCL[0-9]", "^CXCL[0-9]", "^CX3CL", "^XCL"),
  
  # Chemokine receptors
  chemokine_receptors = c("^CCR[0-9]", "^CXCR[0-9]", "^CX3CR"),
  
  # Immunoglobulins and antibodies
  antibodies = c("^IGH", "^IGK", "^IGL", "^IGHV", "^IGKV", "^IGLV"),
  
  # T cell markers and receptors
  tcell = c("^CD3[A-Z]?$", "^CD3[DEG]", "^CD4$", "^CD8[AB]?$", 
            "^CD2$", "^CD5$", "^CD7$", "^CD27$", "^CD28$", "^CD69$", 
            "^CD95$", "^CD122$", "^CD127$", "^CTLA4$", "^PDCD1$", 
            "^LAG3$", "^HAVCR2$", "TIGIT$", "^TCR", "^TRA$", "^TRB$", 
            "^TRG$", "^TRD$"),
  
  # B cell markers
  bcell = c("^CD19$", "^CD20$", "^MS4A1$", "^CD22$", "^PAX5$", "^BCL6$"),
  
  # NK cell markers
  nk = c("^KLRD1$", "^KLRF1$", "^KLRC", "^NCR[0-9]", "^NKG2"),
  
  # Myeloid markers
  myeloid = c("^CD14$", "^CD68$", "^CD163$", "^ITGAM$", "^ITGAX$", 
              "^CSF1R$", "^FCGR3A$", "^CD16$"),
  
  # Complement system
  complement = c("^C1QA$", "^C1QB$", "^C1QC$", "^C1R$", "^C1S$", "^C3$", 
                 "^C4A$", "^C4B$", "^C5$", "^C6$", "^C7$", "^C8", "^C9$", 
                 "^CFH$", "^CFI$", "^CFB$", "^CFD$"),
  
  # MHC/HLA genes
  mhc = c("^HLA-", "^H2-", "^B2M$", "^TAP1$", "^TAP2$", "^PSMB8$", "^PSMB9$"),
  
  # Toll-like receptors
  tlr = c("^TLR[0-9]", "^MYD88$", "^IRAK", "^TRAF"),
  
  # Transcription factors and signaling
  signaling = c("^STAT[0-9]", "^JAK[0-9]", "^SOCS[0-9]", "^IRF[0-9]", 
                "^NFKB", "^REL$", "^RELA$", "^RELB$"),
  
  # T cell subsets
  tsubsets = c("^FOXP3$", "^GATA3$", "^TBX21$", "^RORC$", "^BCL6$"),
  
  # Fc receptors
  fc_receptors = c("^FCGR", "^FCER", "^FCRL"),
  
  # Cytotoxicity
  cytotoxic = c("^PRF1$", "^GZMA$", "^GZMB$", "^GZMH$", "^GNLY$", "^GZMK$"),
  
  # Checkpoint molecules
  checkpoint = c("^PDCD1$", "^CD274$", "^PDCD1LG2$", "^CTLA4$", 
                 "^LAG3$", "^HAVCR2$", "^TIGIT$", "^BTLA$", "^VSIR$"),
  
  # Inflammasome/inflammatory
  inflammasome = c("^NLRP", "^NLRC", "^AIM2$", "^PYCARD$", "^CASP1$", 
                   "^IL1B$", "^IL18$"),
  
  # Adhesion molecules
  adhesion = c("^ICAM", "^VCAM", "^SELL$", "^SELP$", "^SELE$", "^ITGA", "^ITGB"),
  
  # Costimulatory molecules
  costimulatory = c("^CD80$", "^CD86$", "^CD28$", "^ICOS", "^CD40$", "^CD27$")
)

# Create combined immune gene pattern
immune_pattern <- paste(unlist(immune_keywords), collapse = "|")

# ============================================================================
# Set up directories
# ============================================================================
projectdirectory      <- "/data/rdandekar/rprojects/TBM_RNASeq_Analysis"
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
res_df_sig_zi_full$is_immune <- grepl(immune_pattern, 
                                      res_df_sig_zi_full$gene_name, 
                                      ignore.case = FALSE)

# Subset to immune genes only
res_df_immune <- res_df_sig_zi_full[res_df_sig_zi_full$is_immune, ]

# Print summary statistics
cat("\n=== Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(res_df_sig_zi_full), "\n")
cat("Total immune genes identified:", nrow(res_df_immune), "\n")
cat("Significant immune genes (padj < 0.05):", 
    sum(res_df_immune$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated immune genes (padj < 0.05, log2FC > 1):", 
    sum(res_df_immune$padj < 0.05 & res_df_immune$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Downregulated immune genes (padj < 0.05, log2FC < -1):", 
    sum(res_df_immune$padj < 0.05 & res_df_immune$log2FoldChange < -1, na.rm = TRUE), "\n\n")

# Save immune gene list
write.csv(res_df_immune, 
          file = file.path(output_dir, "immune_genes_dge_results.csv"),
          row.names = FALSE, quote = FALSE)

# ============================================================================
# Create Volcano Plots
# ============================================================================
# Define genes to label
genes_to_label <- c('IL1A', 'C1S', 'STAT4', 'IL7', 'IL24', 'TRAF1', 'TRAF2', 
                    'IL3','IL23R', 'TNF', 'IRF6', 'CD27', 'TNFRSF1B', 'IL6ST', 
                    'C4A', 'C4B', 'CXCL1', 'NLRP12', 'CARD', 'TNFRSF13B', 'UBE2N', 'LCK', 'TCR', 'LAT',
                    'FYN', 'IL12RB1', 'BCR', 'CD81', 'PLCG2', 'CARD10', 'STING1', 'TYK2')

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
  title = 'Differentially Expressed Genes: TBM Survivors',
  subtitle = 'p-values capped at -log10(p)=5.0, log2FC capped at ±10 for clarity',
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
  labFace = 'bold',
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
top.table.hclus3.full$is_immune <- grepl(immune_pattern, 
                                         top.table.hclus3.full$gene, 
                                         ignore.case = FALSE)

# Subset to immune genes only
cluster3_immune <- top.table.hclus3.full[top.table.hclus3.full$is_immune, ]

# Print summary statistics
cat("\n=== Immune Gene Summary ===\n")
cat("Total genes in DGE results:", nrow(top.table.hclus3.full), "\n")
cat("Total immune genes identified:", nrow(cluster3_immune), "\n")
cat("Significant immune genes (padj < 0.05):", 
    sum(cluster3_immune$padj < 0.05, na.rm = TRUE), "\n")

sig_immune_only <- cluster3_immune$gene[
  cluster3_immune$padj < 0.05 & !is.na(cluster3_immune$padj)
]

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

cluster3_genes_label <- c('ITGAX', 'CCL5', 'HLA-A', 'CSF2RA', 'TNFSF10', 'BCL6', 
                          'CD14', 'CD40', 'CTLA4', 'IL27RA', 'CCR1', 'CCR6', 
                          'FCGR1A', 'ITGA4', 'STAT5B', 'IL17RA', 'IL1B', 'REL', 
                          "B2M", 'TAP1', 'C1QA', 'C1QC', 'C7', 'CCL16', 'CCL3', 
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
top.table.hclust3.survival$is_immune <- grepl(immune_pattern, 
                                        top.table.hclust3.survival$gene, 
                                         ignore.case = FALSE)

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

sig_immune_only <- cluster3surv_immune$gene[
  cluster3surv_immune$padj < 0.05 & !is.na(cluster3surv_immune$padj)
]

plot_volcano_hc3_surv <- EnhancedVolcano(
  toptable = top.table.hclust3.survival.capped,
  lab = top.table.hclust3.survival.capped$gene,
  x = 'log2FoldChange_capped',
  y = 'padj_capped',
  ylim = c(0, 5),  # Limit y-axis to 0-10
  title = 'Hypoinflammatory Cluster: TBM Survivors',
  subtitle = 'Differentially expressed genes in comparison to cluster nonsurvivors ',
  selectLab = sig_immune_only,
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

