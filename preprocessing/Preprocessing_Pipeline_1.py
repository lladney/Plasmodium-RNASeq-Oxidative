# Plasmodium RNA-seq Analysis for Oxidative Stress Profiling
# Step 1: GEO-to-SRA conversion, metadata extraction, FASTQ download,...
# ...preprocessing (trimming), and quality control (FastQC + MultiQC)

# LIBRARIES
import os                                                       # os = to interact with operating system
from pysradb.sraweb import SRAweb                               # SRAweb = to query metadata from GEO and SRA
import pandas as pd                                             # pandas = to conduct data analysis/manipulation
import subprocess                                               # subprocess = to run external commands from Python (i.e., prefetch, cutadapt)
import shutil                                                   # shutil = to verify tools before running commands

# PATH SETUP FOR CONDA TOOLS
os.environ['PATH']                                              # PATH stores list of directories; os.environ yields current list of folders Python searches
    += (                                                        # Appends additional directories
    f":{os.environ['CONDA_PREFIX']}/bin"                        # When external tools called, Python knows where to find them
                                                                # Component breakdoown:
                                                                # os.environ['CONDA_PREFIX'] = root directory of current Conda environment
                                                                # /bin = targets folder where Conda-installed command-line tools are
                                                                # f = tell Python to treat string as executable
                                                                # Note: need ":" because it separates directories when multiple
    ":/Users/laraladney/Downloads/FastQC"                       # User has to modify this to location where FastQC is installed
)

# CONNECT TO SRA
db = SRAweb()                                                   # Connect to NCBI SRA; assign to db

geo_id = "GSE189034"                                            # Set GEO Series Accession to variable geo_id

try:
    # Try to fetch SRA project(s) from GEO
    gse_info = db.gse_to_srp(geo_id)

    if gse_info is not None and not gse_info.empty:
        print("Returned columns:", gse_info.columns.tolist())
        print(gse_info)

        if 'study_accession' in gse_info.columns:
            srp_id = gse_info['study_accession'].values[0]
            print("Selected SRA Project:", srp_id)

            # Fetch metadata
            metadata = db.sra_metadata(srp_id, detailed=True)
            metadata.to_csv(f"{srp_id}_metadata.csv", index=False)
            print(f"Metadata saved to {srp_id}_metadata.csv")

            # ------------------------------
            # Step 2: Download & Preprocess
            # ------------------------------

            # Check for required tools
            required_tools = ["prefetch", "fasterq-dump", "fastqc", "multiqc", "cutadapt"]
            #
            for tool in required_tools:
                print(f"{tool} path:", shutil.which(tool))

            #
            for tool in required_tools:
                if shutil.which(tool) is None:
                    raise EnvironmentError("Required tool '{tool}' not found in PATH. Please install it.")

            # Extract run accessions (SRR IDs)
            # sra_runs = metadata['run_accession'].tolist()
            sra_runs = [
                "SRR16966869", "SRR16966870", "SRR16966871",
                "SRR16966872", "SRR16966873", "SRR16966874"
            ]

            
            # Create directories
            os.makedirs("raw_data", exist_ok=True)
            os.makedirs("trimmed_data", exist_ok=True)
            os.chdir("raw_data")

            # Download and convert SRA to FASTQ
            for run in sra_runs:
                print(f"Downloading {run}...")
                subprocess.run(["prefetch", run], check=True)
                subprocess.run(["fasterq-dump", run, "--split-files", "--threads", "4"], check=True)

            # Run FastQC on raw FASTQ files
            fastq_files = [f for f in os.listdir(".") if f.endswith(".fastq")]
            print("Running FastQC on raw FASTQ files...")
            subprocess.run(["fastqc"] + fastq_files, check=True)
            
            # Run MultiQC to summarize FastQC reports
            print("Running MultiQC...")
            subprocess.run(["multiqc", "."], check=True)

            # Step 3: Adapter trimming using cutadapt
            print("Trimming adapters using cutadapt...")
            for file in fastq_files:
                trimmed_file = os.path.join("../trimmed_data", f"trimmed_{file}")
                subprocess.run([
                    "cutadapt",
                    "-q", "20",
                    "-m", "30",
                    "-a", "AGATCGGAAGAGC",
                    "-o", trimmed_file,
                    file
                ], check=True)

            print("Preprocessing and QC complete. Trimmed files saved in 'trimmed_data/' directory.")

        else:
            print("Expected column 'study_accession' not found. Available columns printed above.")
    else:
        print(f"No SRA Project found for GEO accession {geo_id}.")

except Exception as e:
    print(f"An error occurred: {e}")


