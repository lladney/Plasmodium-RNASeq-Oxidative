# Plasmodium-RNASeq-Oxidative
A modular pipeline for RNA-seq analysis of oxidative stress responses in Plasmodium, featuring SRA data handling, Cutadapt trimming, Salmon quantification, and DESeq2 differential analysis.

## Summary
1. Preprocessing of raw reads
2. Transcript-level quantification using Salmon
3. Differential gene expression (DGE) analysis using DESeq2 in R

## Project Structure
```
plasmodium-rnaseq-pipeline/
├── preprocessing/                    # Step 1: Download, trim, and QC raw FASTQ files
│   ├── Preprocessing_Pipeline_1.py   # Python script for SRA download + trimming + QC
│   ├── SRP346589_metadata.csv        # Metadata file from SRA/GEO
│   ├── raw_data/                     # Raw FASTQ files + FastQC + MultiQC (gitignored)
│   └── trimmed_data/                 # Trimmed FASTQ reads for quantification
│
├── quantification/                   # Step 2: Transcript quantification using Salmon
│   ├── Alignment_Quantification_2.py # Python script for Salmon indexing & quantification
│   ├── count_matrix.csv              # Gene-level count matrix (for DESeq2)
│   └── metadata.csv                  # Sample metadata (condition/grouping info)
│
├── dgeanalysis/                      # Step 3: Differential gene expression (DGE) analysis in R
│   ├── DGE_Analysis_DESeq2.R         # R script for DESeq2 analysis and visualization
│   ├── deseq2_results.csv            # Full table of DE results
│   ├── PCA_plot.png                  # Principal component analysis of samples
│   ├── heatmap_top20_DEGs.png        # Heatmap of top 20 most significant DEGs
│   └── volcano_plot.png              # Volcano plot highlighting significant DEGs
│
├── reference/                        # Reference files for Salmon quantification
│   ├── transcripts.fa                # Reference transcriptome (not tracked by git)
│   └── tx2gene.csv                   # Transcript-to-gene mapping file
│
├── .gitignore                        # Prevents tracking of large FASTQ, zip, and HTML files
├── requirements.txt                  # Python dependencies (e.g., pysradb, pandas, multiqc)
└── README.md                         # Project overview and usage instructions
```
## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/plasmodium-rnaseq-pipeline.git
cd plasmodium-rnaseq-pipeline
```  

2. Create a Conda environment:
```bash
conda create -n plasmodium_env python=3.10
conda activate plasmodium_env
```

3. Install Python dependencies:
```
pip install -r requirements.txt
```

4. Install the required command-line tools via [bioconda](https://bioconda.github.io/):
   - `fastqc`: Quality control for raw reads
   - `cutadapt`: Adapter trimming
   - `salmon`: Transcript-level quantification
   - `multiqc`: Aggregated QC reporting

   You can install them all at once:
   ```bash
   conda install -c bioconda fastqc cutadapt salmon multiqc

## Running the Pipeline

### Step 1:  *PREPROCESSING*
Go to the preprocessing/ directory and run: 
```python 
Preprocessing_Pipeline_1.py
```
This will perform: 
- GEO/SRA Metadata Fetching
- FASTQ Download
- Adapter Trimming with Cutadapt
- Quality Control with FastQC
- MultiQC Summary

### Step 2:  *QUANTIFICATION*
Go to the quantification/ directory and run:
```python
Alignment_Quantification_2.py
```
This will: 
- Build a Salmon Index
- Quantify transcript expression
- Generate gene-level count matrix and metadata CSV

### Step 3:  *DIFFERENTIAL EXPRESSION ANALYSIS*
Go to the dgeanalysis/ directory and run:
```Rscript
DGE_Analysis_DESeq2.R
```	
This will produce: 
- deseq2_results.csv: Differential expression output
- volcano_plot.png: Volcano plot of DEGs
- heatmap_top20_DEGs.png: Heatmap of top 20 DEGs
- PCA_plot.png: Principal component analysis of samples

## Notes
* Raw data files and large intermediate results are .gitignored
* This pipeline was developed and tested on macOS 10.15 with Conda and R 4.3+
