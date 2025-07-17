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

metadata$condition <- factor(metadata$condition,                 # Condition column in metadata treated as factor (categorical variable)
                             levels = c("control", "dozi_ko"))   # Factor levels: control (ensure it's baseline/reference), dozi_ko (condition compared to control)
                                                                 # Note: DESeq2 will calculate log2 fold change as dozi_ko vs control

stopifnot(all(colnames(counts) == rownames(metadata)))           # Make sure samples names in columns of count matrix match samples names in rows of metadata

# CREATE DESEQ2 DATASET
dds <- DESeqDataSetFromMatrix(countData = round(counts),         # Create DESeq2 dataset object dds; countData = gene count matrix with integer values (round any decimals)
                              colData = metadata,                # Sample metadata; tell DESeq2 if sample belongs to control or dozi_ko
                              design = ~ condition)              # Tell DESeq2 to fit a model where gene expression is explained by condition
                                                                 # Note: Is gene expression significantly different between control and dozi_ko?
                                                                 # dds contains:
                                                                 # 1) Count matrix
                                                                 # 2) Sample metadata
                                                                 # 3) Design 
# RUN DESEQ2
dds <- DESeq(dds)                                                # Run DESeq2 on dds object
                                                                 # DESeq2 functions:
                                                                 # 1) Estimate size factors for normalization
                                                                 # 2) Estimate dispersion (gene expression variance across replicates)
                                                                 # 3) Fit negative binomial generalized linear model (extra variability for noisy count data)

# SHRINK LOG2 FOLD CHANGES 
res <- lfcShrink(dds,                                            # dds with DE analysis
                 coef = "condition_dozi_ko_vs_control",          # Specify which comparison to shrink
                 type = "apeglm")                                # Apply shrinkage 
                                                                 # Note: Shrinkage helps with extreme fold changes, helps with better gene ranking, better visualization, exc.

write.csv(as.data.frame(res), file = output_path)                # Convert DESeq2 results DataFrame to CSV file
                                                                 # CSV contains:
                                                                 # 1) Gene names (row names)
                                                                 # 2) log2FoldChange
                                                                 # 3) pvalue
                                                                 # 4) padj (FDR)
                                                                 # 5) Standard errors

# VOLCANO PLOT (PNG)
png("/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/volcano_plot.png", width = 1200, height = 1000) # Write volcano plot to PNG
EnhancedVolcano(res,                                             # res = contains log2FoldChange, padj
    lab = rownames(res),                                         # lab = gene labels, from rows of res
    x = "log2FoldChange",                                        # x-axis (magnitude of expression change)
    y = "padj",                                                  # y-axis (adjusted p-value/significance)
    pCutoff = 0.05,                                              # Horizontal threshold (adjP/FDR = 0.05)
    FCcutoff = 1.0,                                              # Vertical threshold (log2FC +/- 1 (fold change of 2))
    title = "DOZI KO vs Control (40h Stress)")                   # Title
dev.off()                                                        # Close PNG and write plot

# HEATMAP OF TOP 20 DEGS (PNG) 
top <- head(order(res$padj,                                      # Rank genes by padj, smallest (most significant) to largest
                  na.last = NA),                                 # Remove NA values (missing)
            20)                                                  # Get 20 most significant genes
vsd <- vst(dds,                                                  # Apply variance-stabilizing transformation to raw counts from dds
           blind = FALSE)                                        # Use experimental design (control vs dozi_ko)

png("/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/heatmap_top20_DEGs.png", width = 1000, height = 800) # Write heatmap to PNG
pheatmap(assay(vsd)[top, ],                                      # Get VST-transformed expression matrix for top 20 genes
         cluster_rows = TRUE,                                    # Cluster genes based on expression profiles
         show_rownames = TRUE,                                   # Show gene names on heatmap rows
         cluster_cols = TRUE,                                    # Cluster samples by expression similarity
         annotation_col = metadata)                              # Use sample metadata to annotate columns by group
dev.off()                                                        # Close PNG and write plot

# PCA PLOT (PNG)
png("/Users/laraladney/Documents/plasmodium-rnaseq-pipeline/dgeanalysis/PCA_plot.png", width = 1000, height = 800) # Write PCA plot to PNG
plotPCA(vsd, intgroup = "condition")                             # Perform PCA on VST-transformed expression matrix, coloring samples by experimental condition
dev.off()                                                        # Close PNG and write plot
