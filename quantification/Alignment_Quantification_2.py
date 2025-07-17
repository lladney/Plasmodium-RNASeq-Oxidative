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
    r1 = fastq_dir / f"trimmed_{sample}_1.fastq"                  # Full path to read 1 .fastq
    r2 = fastq_dir / f"trimmed_{sample}_2.fastq"                  # Full path to read 2 .fastq
    output_dir = salmon_out_dir / sample                          # Create output directory for sample's Salmon results

    if not r1.exists() or not r2.exists():                        # Check that both .fastq files exist
        raise FileNotFoundError(f"Could not find paired-end FASTQ files: {r1}, {r2}") # Notify user if not

    quant_cmd = [                                                 # Salmon quantification command
        "salmon", "quant",                                        # Run salmon quant subcommand
        "-i", str(index_dir),                                     # Input: path to Salmon transcriptome index
        "-l", "A",                                                # Library type: A, for automatic detection of paired-end reads
        "-1", str(r1),                                            # Path to read 1 .fastq (forward reads)
        "-2", str(r2),                                            # Path to read 2 .fastq (reverse reads)
        "-p", "8",                                                # Use 8 threads for parallel processing
        "--validateMappings",                                     # Verify that mappings are consistent with transcript sequences
        "-o", str(output_dir)                                     # Output directory
    ]

    subprocess.run(quant_cmd, check=True)                         # Execute Salmon command; stop if command fails
    print(f"Quantification complete for {sample}")                # Notify user

# BUILD METADATA

print("Generating metadata.csv...")                               # Notify user

metadata = pd.DataFrame({                                         # Create pandas DataFrame
    "condition": [sample_conditions[sample]                       # Loops through each sample in samples_names and looks up condition in sample_conditions
                  for sample in sample_names]
}, index=sample_names)                                            # Set row names (index) of DataFrame to match sample names
metadata.index.name = "sample"                                    # Set index column name to "sample"
metadata_path = quant_dir / "metadata.csv"                        # Create path where CSV will save
metadata.to_csv(metadata_path)                                    # Save metadata DataFrame as CSV in created path

print(f"Metadata saved to {metadata_path}")                       # Notify user

# BUILD COUNT MATRIX USING TX2GENE

print("Building gene-level count matrix using tx2gene.csv...")    # Notify user
                                                                  # Note: DESeq2 runs at gene level -> need to aggregate transcript counts by gene
tx2gene_df = pd.read_csv(tx2gene_path,                            # Read CSV into a pandas DataFrame
                         names=["TXNAME", "GENEID"])              # Assign column names

count_dict = {}                                                   # Create dictionary to store gene-level counts for each sample
gene_set = set()                                                  # Initialize empty set to collect unique gene IDs for each sample

for sample in sample_names:                                       # Loop through each sample in samples_names
    quant_file = salmon_out_dir / sample / "quant.sf"             # Create full path to Salmon output file
    quant_df = pd.read_csv(quant_file,                            # Load quant.sf into pandas DataFrame, quant_df
                                                                  # Note: quant.sf = tab-delimited
                           sep="\t")[["Name", "NumReads"]]        # Keep only these two columns

    merged_df = quant_df.merge(tx2gene_df,                        # Merge quant_df (only transcript IDs, read count) with tx2gene_df (maps transcript IDs to gene IDs)
                               left_on="Name",                    # To get transcript IDs in quant_df, look at "Name" column
                               right_on="TXNAME")                 # To get transcript IDs in tx2gene_df, look at "TXNAME" column

    gene_counts = merged_df.groupby("GENEID")["NumReads"].sum()   # Group by GENEID and add all transcripts that map to that same GENEID (multiple transcripts per gene)
                                                                  # gene_counts: index = gene ID, value = summed read count
    count_dict[sample] = gene_counts                              # Store gene-level counts from gene_counts in count_dict dictionary under sample name key
    gene_set.update(gene_counts.index)                            # Add all gene IDs in sample quantification to gene_set  
                                                                  # Note: .update avoids duplicates; gene_set = all unique gene IDs across all samples

all_genes = sorted(gene_set)                                      # Sort the full list of gene IDs from the set alphabetically
count_matrix = pd.DataFrame(index=all_genes)                      # Initialize new empty DataFrame where the row index is a sorted list of gene IDs

for sample in sample_names:                                       # Loop through each sample in sample_names
    count_matrix[sample] = count_dict[sample].reindex(all_genes, fill_value=0) # Adds column to count_matrix with gene counts from count_dict
                                                                  # Component Breakdown:
                                                                  # count_dict[sample] = panda Series of gene-level counts for one sample (index = gene IDs, only genes with reads in that sample are present)
                                                                  # reindex(all_genes, fill_value=0) = reorder and align panda Series to full list of genes; missing genes (not in this sample) get a value of 0
                                                                  # Note: ensures that every gene appears in every sample even if count is zero
                                                                  # count_matrix[sample] = assign reindexed Series as a new column in count_matrix labeled with sample name

count_matrix.index.name = "gene"                                  # Set name of index column in DataFrame to "gene"
count_matrix_path = quant_dir / "count_matrix.csv"                # Create path to where matrix will save
count_matrix.to_csv(count_matrix_path)                            # Convert DataFrame to CSV
print(f"Gene-level count matrix saved to {count_matrix_path}")    # Notify user
