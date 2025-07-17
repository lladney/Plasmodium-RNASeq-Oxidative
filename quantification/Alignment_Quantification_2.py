# Plasmodium RNA-seq Analysis for Oxidative Stress Profiling
# Step 2: Quantify with Salmon and summarize results for DE analysis

# LIBRARIES
import os                                                         # os = to interact with operating system
import pandas as pd                                               # pandas = to conduct data analysis/manipulation
import subprocess                                                 # subprocess = to run external commands from Python (i.e., Salmon)
import shutil                                                     # shutil = to verify tools before running commands
from pathlib import Path                                          # Path = to get file locations, create directories, exc.

# ------------ CONFIGURATION ------------ #

# Set working directories
root_dir = Path(__file__).resolve().parents[1]  # Project root
quant_dir = root_dir / "quantification"
fastq_dir = root_dir / "preprocessing" / "trimmed_data" # Adjust as needed
index_dir = quant_dir / "salmon_index"
salmon_out_dir = quant_dir / "salmon_output"
salmon_out_dir.mkdir(parents=True, exist_ok=True)

# Input files and sample names
sra_runs = [
    "SRR16966869", "SRR16966870", "SRR16966871",
    "SRR16966872", "SRR16966873", "SRR16966874"
]
sample_names = sra_runs

# Mapping sample names to experimental conditions
sample_conditions = {
    "SRR16966869": "control",
    "SRR16966870": "control",
    "SRR16966871": "control",
    "SRR16966872": "dozi_ko",
    "SRR16966873": "dozi_ko",
    "SRR16966874": "dozi_ko"
}

# Reference transcriptome
transcriptome_fasta = root_dir / "reference" / "Pfalciparum_transcripts.fa"

# ------------ BUILD TX2GENE ------------ #

tx2gene_path = root_dir / "reference" / "tx2gene.csv"
if not tx2gene_path.exists():
    print("📄 Generating tx2gene.csv from FASTA headers...")
    tx2gene = []
    with open(transcriptome_fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                header = line.strip().split()[0][1:]
                gene_id = header.split(".")[0]
                tx2gene.append((header, gene_id))
    pd.DataFrame(tx2gene, columns=["TXNAME", "GENEID"]).to_csv(tx2gene_path, index=False)
    print(f"Saved: {tx2gene_path}")

# ------------ RUN SALMON QUANTIFICATION ------------ #

print("🧬 Running Salmon quantification...")
for sample in sample_names:
    r1 = fastq_dir / f"trimmed_{sample}_1.fastq"
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

print("📄 Generating metadata.csv...")
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
