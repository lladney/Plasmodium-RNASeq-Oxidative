# Plasmodium RNA-seq Analysis for Oxidative Stress Profiling
# Step 2: Quantify with Salmon and summarize results for DE analysis

# LIBRARIES
import os                                                         # os = to interact with operating system
import pandas as pd                                               # pandas = to conduct data analysis/manipulation
import subprocess                                                 # subprocess = to run external commands from Python (i.e., Salmon)
import shutil                                                     # shutil = to verify tools before running commands
from pathlib import Path                                          # Path = to get file locations, create directories, exc.

# SET WORKING DIRECTORIES
root_dir = Path(__file__).resolve().parents[1]                    # Get absolute path to root directory
                                                                  # __file__ = path to current script
                                                                  # resolve() = full absolute path
                                                                  # .parents[1] = go up one level from project folder
quant_dir = root_dir / "quantification"                           # Define folder for quantification
fastq_dir = root_dir / "preprocessing" / "trimmed_data"           # Define folder where trimmed .fastq files saved (in script 01)
index_dir = quant_dir / "salmon_index"                            # Define folder for Salmon index in quantification folder
salmon_out_dir = quant_dir / "salmon_output"                      # Define folder for Salmon output in quantification folder
salmon_out_dir.mkdir(parents=True, exist_ok=True)                 # Create salmon_outut and parents if doesn't exist yet

# INPUT FILES AND SAMPLE NAMES
sra_runs = [                                                      # List of SRR IDs (sra_runs)
    "SRR16966869", "SRR16966870", "SRR16966871",
    "SRR16966872", "SRR16966873", "SRR16966874"
]
sample_names = sra_runs                                           # Make second list but assign nane sample_names
                                                                  # Preserve raw sequencing run IDs (sra_runs) and sample_names if renamed
# MAP SAMPLE NAMES TO EXPERIMENTAL CONDITIONS 
sample_conditions = {                                             # SRR ID either maps to "control" or "dozi_ko"
    "SRR16966869": "control",    
    "SRR16966870": "control",
    "SRR16966871": "control",
    "SRR16966872": "dozi_ko",
    "SRR16966873": "dozi_ko",
    "SRR16966874": "dozi_ko"
}

# REFERENCE TRANSCRIPTOME
transcriptome_fasta = root_dir / "reference" / "Pfalciparum_transcripts.fa" # Path to reference transcriptome fasta

# BUILD TX2GENE
tx2gene_path = root_dir / "reference" / "tx2gene.csv"             # Set path to transcript-to-gene mapping
if not tx2gene_path.exists():                                     # Don't make file if already exists 
    print("Generating tx2gene.csv from FASTA headers...")         # Notify user
    tx2gene = []                                                  # Create empty list                       
    with open(transcriptome_fasta, 'r') as fasta:                 # Open FASTA containing transcriptome sequences
        for line in fasta:                                        # Loop through each line in .fa
            if line.startswith(">"):                              # Just process header lines
                header = line.strip().split()[0][1:]              # Get transcript ID (isoform) and remove carat
                gene_id = header.split(".")[0]                    # Get gene ID (gene root) and remove the isoform value
                tx2gene.append((header, gene_id))                 # Add the transcript and gene IDs to the list
    pd.DataFrame(tx2gene,                                         # tx2gene is list of tuples
                 columns=["TXNAME", "GENEID"]).to_csv(tx2gene_path, # Name columns for transcript and gene IDs and write to CSV
                                            index=False)          # Don't write row numbers
    print(f"Saved: {tx2gene_path}")                               # Notify user

# RUN SALMON QUANTIFICATION
print("Running Salmon quantification...")                         # Notify user
for sample in sample_names:                                       # Loop through list of sample IDs
    r1 = fastq_dir / f"trimmed_{sample}_1.fastq"                  # 
    r2 = fastq_dir / f"trimmed_{sample}_2.fastq"
    output_dir = salmon_out_dir / sample

    if not r1.exists() or not r2.exists():
        raise FileNotFoundError(f"Could not find paired-end FASTQ files: {r1}, {r2}")

    quant_cmd = [
        "salmon", "quant",
        "-i", str(index_dir),
        "-l", "A",
        "-1", str(r1),
        "-2", str(r2),
        "-p", "8",
        "--validateMappings",
        "-o", str(output_dir)
    ]

    subprocess.run(quant_cmd, check=True)
    print(f"Quantification complete for {sample}")

# ------------ BUILD METADATA ------------ #

print("Generating metadata.csv...")
# Build and save metadata with sample as index (for DESeq2)
# Build and save metadata with sample as index (for DESeq2)
metadata = pd.DataFrame({
    "condition": [sample_conditions[sample] for sample in sample_names]
}, index=sample_names)
metadata.index.name = "sample"
metadata_path = quant_dir / "metadata.csv"
metadata.to_csv(metadata_path)

print(f"Metadata saved to {metadata_path}")

# ------------ BUILD COUNT MATRIX (using tx2gene) ------------ #

print("Building gene-level count matrix using tx2gene.csv...")
tx2gene_df = pd.read_csv(tx2gene_path, names=["TXNAME", "GENEID"])

count_dict = {}
gene_set = set()

for sample in sample_names:
    quant_file = salmon_out_dir / sample / "quant.sf"
    quant_df = pd.read_csv(quant_file, sep="\t")[["Name", "NumReads"]]

    # Merge with tx2gene mapping
    merged_df = quant_df.merge(tx2gene_df, left_on="Name", right_on="TXNAME")

    # Group by GENEID and sum counts
    gene_counts = merged_df.groupby("GENEID")["NumReads"].sum()
    count_dict[sample] = gene_counts
    gene_set.update(gene_counts.index)

# Combine all samples into one DataFrame
all_genes = sorted(gene_set)
count_matrix = pd.DataFrame(index=all_genes)

for sample in sample_names:
    count_matrix[sample] = count_dict[sample].reindex(all_genes, fill_value=0)

count_matrix.index.name = "gene"
count_matrix_path = quant_dir / "count_matrix.csv"
count_matrix.to_csv(count_matrix_path)
print(f"Gene-level count matrix saved to {count_matrix_path}")
