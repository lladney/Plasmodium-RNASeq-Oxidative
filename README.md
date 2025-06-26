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

## Dataset

This project uses bulk RNA-seq data from [SRP346589](https://www.ncbi.nlm.nih.gov/sra?term=SRP346589), profiling *Plasmodium falciparum* parasites exposed to oxidative stress.

- **Title:** Transcriptional response of *Plasmodium falciparum* to peroxide antimalarials  
- **SRA Accession:** SRP346589  
- **Organism:** *Plasmodium falciparum* (3D7 strain)  
- **Conditions:** DMSO (control) vs. peroxide antimalarials (e.g., OZ277, OZ439) at multiple timepoints
- **SRA Runs Used:**  
  - **SRR16966869** – DMSO control, 3 hr  
  - **SRR16966870** – OZ277 treatment, 3 hr  
  - **SRR16966871** – OZ439 treatment, 3 hr  
  - **SRR16966872** – DMSO control, 6 hr  
  - **SRR16966873** – OZ277 treatment, 6 hr  
  - **SRR16966874** – OZ439 treatment, 6 hr    
- **Rationale:** Chosen to explore oxidative stress–induced gene expression, the dataset includes high-quality paired-end bulk RNA-seq suitable for differential expression and functional enrichment analysis. Sample metadata is organized in `SRP346589_metadata.csv`.

### Citation

**PubMed (APA format):**  
Bruning-Richardson, A., Coomes, D., Crespo, M. P., Bottrill, A. R., Wilkinson, S. R., & Ward, S. A. (2022). *The transcriptional response of Plasmodium falciparum to peroxide antimalarials*. *Frontiers in Cellular and Infection Microbiology, 12*, 841957. https://doi.org/10.3389/fcimb.2022.841957

**BibTeX:**
```bibtex
@article{BruningRichardson2022Peroxide,
  author  = {Bruning-Richardson, A. and Coomes, D. and Crespo, M. P. and Bottrill, A. R. and Wilkinson, S. R. and Ward, S. A.},
  title   = {The transcriptional response of Plasmodium falciparum to peroxide antimalarials},
  journal = {Frontiers in Cellular and Infection Microbiology},
  year    = {2022},
  volume  = {12},
  pages   = {841957},
  doi     = {10.3389/fcimb.2022.841957}
}
```
## Notes
* Raw data files and large intermediate results are excluded via ```.gitignore```
* This pipeline was developed and tested on macOS 10.15 with Conda and R 4.3+
