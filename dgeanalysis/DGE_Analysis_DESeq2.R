# Set default CRAN mirror
cran_mirror <- "https://cloud.r-project.org"

# Ensure mvtnorm is available (CRAN package)
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    install.packages("mvtnorm", repos = cran_mirror)
}

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_mirror)
}

# Install required Bioconductor packages if missing
bioc_pkgs <- c("DESeq2", "apeglm", "pheatmap", "EnhancedVolcano")
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg)
    }
}

# Load libraries
library(DESeq2)
library(apeglm)
library(pheatmap)
library(EnhancedVolcano)

# Set absolute paths
count_path <- "/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/quantification/count_matrix.csv"
meta_path <- "/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/quantification/metadata.csv"
output_path <- "/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/deseq2_results.csv"

# Load data
counts <- read.csv(count_path, row.names = 1)
metadata <- read.csv(meta_path, row.names = 1)

# Ensure condition is a factor and ordered
metadata$condition <- factor(metadata$condition, levels = c("control", "dozi_ko"))

# Ensure sample names match
stopifnot(all(colnames(counts) == rownames(metadata)))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Shrink log2 fold changes
res <- lfcShrink(dds, coef = "condition_dozi_ko_vs_control", type = "apeglm")

# Save results
write.csv(as.data.frame(res), file = output_path)

# Volcano plot (save as PNG)
png("/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/volcano_plot.png", width = 1200, height = 1000)
EnhancedVolcano(res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1.0,
    title = "DOZI KO vs Control (40h Stress)")
dev.off()

# Heatmap of top 20 DEGs
top <- head(order(res$padj, na.last = NA), 20)
vsd <- vst(dds, blind = FALSE)

png("/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/heatmap_top20_DEGs.png", width = 1000, height = 800)
pheatmap(assay(vsd)[top, ],
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         annotation_col = metadata)
dev.off()

# PCA plot (save as PNG)
png("/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/PCA_plot.png", width = 1000, height = 800)
plotPCA(vsd, intgroup = "condition")
dev.off()

