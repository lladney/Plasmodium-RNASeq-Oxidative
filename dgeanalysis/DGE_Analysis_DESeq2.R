# Plasmodium RNA-seq Analysis for Oxidative Stress Profiling
# Step 3: DEG Analysis and Visualization in R

# SET CRAN MIRROR
cran_mirror <- "https://cloud.r-project.org"                      # cran_mirror holds URL for CRAN mirror (server where R packages downloaded from)

# MVTNORM INSTALLED? 
if (!requireNamespace("mvtnorm", quietly = TRUE)) {               # Check if mvtnorm installed (true)
    install.packages("mvtnorm", repos = cran_mirror)              # Install mvtnorm from CRAN mirror if not (false)      
}                                                                 # Note: need mvtnorm for log2 fold change shrinkage in DESeq2

# BIOCMANAGER INSTALLED?
if (!requireNamespace("BiocManager", quietly = TRUE)) {           # Check if BiocManager installed (true)
    install.packages("BiocManager", repos = cran_mirror)          # Install BiocManager from CRAN mirror if not (false)  
}                                                                 # Note: need BiocManager to install/update packages from Bioconductor

# BIOCONDUCTOR PACKAGES INSTALLED?
                                                                  # Create vector of required packages
bioc_pkgs <- c("DESeq2",                                          # DESeq2 = differential expression analysis
               "apeglm",                                          # apeglm = shrinks log2 fold changes
               "pheatmap",                                        # pheatmap = creates heatmaps
               "EnhancedVolcano")                                 # EnhancedVolcano = creates volcano plots
for (pkg in bioc_pkgs) {                                          # Loop through each package in the vector
    if (!requireNamespace(pkg, quietly = TRUE)) {                 # Check if Bioconductor package installed (true)
        BiocManager::install(pkg)                                 # Install package from BiocManager if not (false)   
    }
}

# LIBRARIES
library(DESeq2)                                                   # Load installed packages into memory
library(apeglm)                                                
library(pheatmap)
library(EnhancedVolcano)

# SET ABSOLUTE PATHS (MODIFY AS NEEDED)
count_path <- "/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/quantification/count_matrix.csv" # Path to gene-level count matrix (quantification script)
meta_path <- "/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/quantification/metadata.csv"      # Path to sample metadata (quantification script)
output_path <- "/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/deseq2_results.csv" # Path to where DESeq2 results will be written as CSV

# LOAD DATA
counts <- read.csv(count_path,                                   # Load gene count matrix where: rows = genes, columns = samples
                   row.names = 1)                                # Tell R that first column of gene IDs are row names not data
metadata <- read.csv(meta_path,                                  # Load metadata mapping samples to experimental condition
                     row.names = 1)                              # Tell R that first column of SRR IDs are row names not data

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

