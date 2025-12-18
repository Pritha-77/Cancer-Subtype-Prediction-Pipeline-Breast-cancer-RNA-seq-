setwd("E:/cancer subtype prediction")

library(GEOquery)
library(data.table)
library(edgeR)
library(limma)
library(dplyr)
library(xtable)
library(rmeta)
library(caret)
library(genefu)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(GSVA)
library(pheatmap)
#devtools::install_github("ccchang0111/PAM50")
library(PAM50)
#___________________Data Preparation (GSE209998)__________________#


# Fetch data from GEO
gse_meta <- getGEO("GSE209998")
gse_meta <- pData(gse_meta[[1]])
rownames(gse_meta) <- gse_meta$description.1

# Read count data
gse_count <- fread("GSE209998_AUR_129_raw_counts.txt", data.table = F)

# Verify sample matching
if(!identical(rownames(gse_meta), colnames(gse_count)[-1])) {
  stop("Sample names don't match between metadata and count data")
}

# Process count table
gse_count <- gse_count[!duplicated(gse_count[[1]]), ]
rownames(gse_count) <- gse_count[[1]]
gse_count <- gse_count[, -1]

# Filter low-expressed genes
perc_keep <- 0.8  # Keep genes expressed in 80% of samples
gene_keep <- rowSums(gse_count > 0) >= ceiling(perc_keep * ncol(gse_count))
count_tbl_low_rm <- gse_count[gene_keep, ]

# Create DGE list and normalize
dge_gse <- DGEList(counts=count_tbl_low_rm, samples = gse_meta)
dge_gse <- dge_gse[, dge_gse$samples$tissue.ch1 == "Breast"]
dge_gse <- calcNormFactors(dge_gse, method = "TMM")
dge_gse_v <- voom(dge_gse, plot=TRUE)


#_________________________Preparing Gene Annotations____________________#

# Get ENTREZ IDs
entrezid <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  keys = rownames(dge_gse_v),
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
)

# Clean up annotations
entrezid <- na.omit(entrezid)
entrezid <- entrezid[!duplicated(entrezid[[1]]) & !duplicated(entrezid[[2]]), ]
rownames(entrezid) <- entrezid[[2]]

# Match genes
com_gene <- intersect(entrezid[[2]], rownames(dge_gse_v))
entrezid <- entrezid[com_gene, c(2, 2, 1)]
colnames(entrezid) <- c("probe", "Gene.Symbol", "EntrezGene.ID")
rownames(entrezid) <- entrezid[[1]]

# Update expression data
dge_gse_v_entrez <- dge_gse_v[rownames(entrezid), ]
dge_gse_v_entrez$genes <- entrezid


#_____________________Subtype Prediction Methods___________________#

#Method 1: PAM50 Classification
#The PAM50 classifier from the PAM50 package uses a 50-gene signature to predict breast cancer subtypes.
# Format expression matrix
expr_dge_gse_v_entrez <- as.data.frame(dge_gse_v_entrez$E)
rownames(expr_dge_gse_v_entrez) <- entrezid$EntrezGene.ID

# Predict subtypes

library(conflicted)
conflicted::conflict_prefer("select", "dplyr") # Always use dplyr's select

df.pred <- PAM50(expr_dge_gse_v_entrez, cutoff = 0)

# Save results
fwrite(as.data.table(df.pred), "subtype_pred_pam50.csv")

# Visualize results
pdf("heatmaptype_pred_pam50_expression.pdf", width = 7, height = 9)
plotPAM50(as.data.frame(dge_gse_v_entrez$E), df.pred$PAM50)
dev.off()

#Method 2: genefu Package


# Transpose expression matrix
expr_dge_gse_v_entrez_t <- t(dge_gse_v_entrez$E)
 
# Load PAM50 data
data(pam50.robust)
 
# Predict subtypes
subtype_pred_pam50 <- molecular.subtyping(
    sbt.model = "pam50",
    data = expr_dge_gse_v_entrez_t,
    annot = entrezid,
    do.mapping = TRUE
)
 
# Process results
subtype_pred_pam50_df <- as.data.frame(subtype_pred_pam50$subtype.proba)
subtype_pred_pam50_df$Subtype <- as.character(subtype_pred_pam50$subtype)
 
# Save results
fwrite(subtype_pred_pam50_df, "subtype_pred_pam50_df.csv")
 
# Visualize results
pdf("heatmaptype_pred_pam50.pdf", width = 7, height = 9)
pheatmap(
    subtype_pred_pam50$subtype.proba,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE
)
dev.off()

# Transpose expression matrix
expr_dge_gse_v_entrez_t <- t(dge_gse_v_entrez$E)

# Load PAM50 data
data(pam50.robust)

# Predict subtypes
subtype_pred_pam50 <- molecular.subtyping(
  sbt.model = "pam50",
  data = expr_dge_gse_v_entrez_t,
  annot = entrezid,
  do.mapping = TRUE
)

# Process results
subtype_pred_pam50_df <- as.data.frame(subtype_pred_pam50$subtype.proba)
subtype_pred_pam50_df$Subtype <- as.character(subtype_pred_pam50$subtype)

# Save results
fwrite(subtype_pred_pam50_df, "subtype_pred_pam50_df.csv")

# Visualize results
pdf("heatmaptype_pred_pam50.pdf", width = 7, height = 9)
pheatmap(
  subtype_pred_pam50$subtype.proba,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE
)
dev.off()


#Method 3: GSVA Approach
# Define custom gene signatures
# Method 3: GSVA / ssGSEA Approach

library(GSVA)
library(pheatmap)
library(data.table)

# Define custom gene signatures
# (Replace character() with actual gene vectors)
gene_sig_ls <- list(
  gene_basal1   = c("MRPS14", "DIP2B"),
  gene_basal2   = c("AKR1C1", "AMIGO2"),
  gene_luminal1 = c("ABCA3", "ABCG1"),
  gene_luminal2 = c("ADAMTS5", "ALCAM")
)

# Run ssGSEA for single-sample profiling
ssgseapar <- ssgseaParam(
  exprData = dge_gse_v$E,
  geneSets = gene_sig_ls
)

gsva_pred <- gsva(ssgseapar)

# ---- Alternative GSVA (commented, valid syntax) ----

# Process and save results
gsva_pred_t <- as.data.frame(t(gsva_pred))
fwrite(gsva_pred_t, file = "gsva_pred_t.csv", row.names = TRUE)

# Visualize results
pdf("heatmap_subtype_pred_ssgsea.pdf", width = 7, height = 9)
pheatmap(
  gsva_pred_t,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE
)
dev.off()
