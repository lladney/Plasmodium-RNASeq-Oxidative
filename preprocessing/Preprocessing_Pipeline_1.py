# Plasmodium RNA-seq Analysis for Oxidative Stress Profiling
# Step 1: GEO-to-SRA conversion, metadata extraction, FASTQ download,...
# ...preprocessing (trimming), and quality control (FastQC + MultiQC)

# LIBRARIES
import os                                                        # os = to interact with operating system
from pysradb.sraweb import SRAweb                                # SRAweb = to query metadata from GEO and SRA
import pandas as pd                                              # pandas = to conduct data analysis/manipulation
import subprocess                                                # subprocess = to run external commands from Python (i.e., prefetch, cutadapt)
import shutil                                                    # shutil = to verify tools before running commands

# PATH SETUP FOR CONDA TOOLS
os.environ['PATH']                                               # PATH stores list of directories; os.environ yields current list of folders Python searches
    += (                                                         # Appends additional directories
    f":{os.environ['CONDA_PREFIX']}/bin"                         # When external tools called, Python knows where to find them
                                                                 # Component breakdoown:
                                                                 # os.environ['CONDA_PREFIX'] = root directory of current Conda environment
                                                                 # /bin = targets folder where Conda-installed command-line tools are
                                                                 # f = tell Python to treat string as executable
                                                                 # Note: need ":" because it separates directories when multiple
    ":/Users/laraladney/Downloads/FastQC"                        # User has to modify this to location where FastQC is installed
)

# CONNECT TO SRA
db = SRAweb()                                                    # Connect to NCBI SRA; assign to db

geo_id = "GSE189034"                                             # Set GEO Series Accession to geo_id

try:                                                             # Let's see if this works
    # Try to fetch SRA project(s) from GEO
    gse_info = db.gse_to_srp(geo_id)                             # Query SRA database (db) to convert GEO accession (geo_id) into SRA project accessions (SRP IDs)
                                                                 # Store SRP IDs in gse_info
    
    if gse_info is not None and not gse_info.empty:              # Must pass 2 conditions:
                                                                 # 1) gse_info is not None (no SRP IDs)
                                                                 # 2) gse_info is not empty (bare DataFrame)
        print("Returned columns:", gse_info.columns.tolist())    # List out column names in DataFrame
        print(gse_info)                                          # Prints DataFrame with SRP IDs

        if 'study_accession' in gse_info.columns:                # Check whether study_accession column exists in gse_info DataFrame
                                                                 # If yes, continue with metadata extraction/download
                                                                 # If no, bypass and print error message to user
            srp_id = gse_info['study_accession'].values[0]       # Extract SRP ID from study_accession column and store it srp_id
                                                                 # Component breakdown:
                                                                 # gse_info['study_accession'] = get the whole study_accession column
                                                                 # .values = convert to NumPy array
                                                                 # [0] = access first entry
            print("Selected SRA Project:", srp_id)               # Print SRP ID so user can confirm    

            # GET METADATA
            metadata = db.sra_metadata(srp_id,                   # Get metadata from SRP ID
                                       detailed=True)            # Return metadata with detailed output  
            metadata.to_csv(f"{srp_id}_metadata.csv",            # Save SRA metadata DataFrame to CSV with SRP ID in name
                            index=False)                         # Don't write row numbers
            print(f"Metadata saved to {srp_id}_metadata.csv")    # Inform user

            # DOWNLOAD AND PREPROCESS

            # CHECK FOR REQUIRED TOOLS
            required_tools = ["prefetch",                        # prefetch = downloads .sra files from NCBI SRA
                              "fasterq-dump",                    # fasterq-dump = convert .sra -> .fastq
                              "fastqc",                          # fastqc = quality control on .fastq files
                              "multiqc",                         # multiqc = put QC from fastqc into HTML summary
                              "cutadapt"]                        # cutadapt = trim adapter sequences, low-quality bases, short reads (clean .fastq files)
            
            for tool in required_tools:                          # Loop through each tool in the required list
                print(f"{tool} path:", shutil.which(tool))       # Return path to tool if available

            for tool in required_tools:                          # Loop through each tool in the required list
                if shutil.which(tool) is None:                   # If tool is not available, stop
                    raise EnvironmentError("Required tool '{tool}' not found in PATH. Please install it.") # Returns error to let user know which tool(s) need to be installed

            # EXTRACT RUN ACCESSIONS (SRR IDS) 
            # sra_runs = metadata['run_accession'].tolist()      # Can get SRR IDs this way but want to run pipeline with smaller subset of data for time purposes
            sra_runs = [                                         # Manual list of SRR IDs from chosen SRA
                "SRR16966869", "SRR16966870", "SRR16966871",
                "SRR16966872", "SRR16966873", "SRR16966874"
            ]
            
            # CREATE DIRECTORIES
            os.makedirs("raw_data",                              # Make raw_data directory
                        exist_ok=True)                           # If already exists, move on
            os.makedirs("trimmed_data",                          # Make trimmed_data directory
                        exist_ok=True)                           # If already exists, move on
            os.chdir("raw_data")                                 # Change working directory to raw_data

            # DOWNLOAD AND CONVERT SRA TO FASTQ
            for run in sra_runs:                                 # Loop through each SRR ID in sra_runs list
                print(f"Downloading {run}...")                   # Print to user which SRR ID currently being processed
                subprocess.run(["prefetch", run],                # Run prefetch to download .sra file for that SRR ID
                               check=True)                       # If download fails, stop
                
                subprocess.run(["fasterq-dump", run,             # Convert downloaded .sra -> .fastq 
                                "--split-files",                 # Split paired-end reads into 2 separate files
                                "--threads", "4"],               # Let's speed this up by using 4 threads
                               check=True)                       # Stop if conversion fails

            # RUN FASTQC ON RAW FASTQ FILES 
            fastq_files = [f for f in os.listdir(".")            # Go through each file in current working directory (".") and list it out
                           if f.endswith(".fastq")]              # Just list out the files ending in .fastq
            print("Running FastQC on raw FASTQ files...")        # Notify user
            subprocess.run(["fastqc"] + fastq_files,             # Run fastqc on all .fastq files found
                           check=True)                           # If fastqc fails, stop
            
            # RUN MULTIQC TO SUMMARIZE FASTQC REPORTS 
            print("Running MultiQC...")                          # Notify user
            subprocess.run(["multiqc", "."],                     # Run multiqc on all files in current working directory (we are still in raw_data/) to make HTML summary
                           check=True)                           # If multiqc fails, stop

            # ADAPTER TRIMMING WITH CUTADAPT
            print("Trimming adapters using cutadapt...")         # Notify user
            for file in fastq_files:                             # Loop through each file in fastq_files
                trimmed_file = os.path.join("../trimmed_data",   # Go to trimmed_data/ and make output path for trimmed file
                                            f"trimmed_{file}")
                subprocess.run([                                 # Run cutadapt
                    "cutadapt",
                    "-q", "20",                                  # Trim low-quality reads from both ends (Phred < 20 (1% error rate))
                    "-m", "30",                                  # Post-trimming, remove reads shorter than 30 bases (probably artifacts)
                    "-a", "AGATCGGAAGAGC",                       # Remove this Illumina adapter sequence to avoid adapter contamination
                    "-o", trimmed_file,                          # Save the output file to trimmed_file
                    file                                         # Input .fastq file being trimmed
                ], check=True)                                   # If cutadapt fails, stop

            print("Preprocessing and QC complete. Trimmed files saved in 'trimmed_data/' directory.") # Notify user

        else:                                                    # If gse_info doesn't have a study_accession column, throws error to the user
            print("Expected column 'study_accession' not found. Available columns printed above.")
    else:                                                        # If gse_info = None or no DataFrame, throws error to the user
        print(f"No SRA Project found for GEO accession {geo_id}.")

except Exception as e:                                           # If the try block fails, throw error to the user
    print(f"An error occurred: {e}")
