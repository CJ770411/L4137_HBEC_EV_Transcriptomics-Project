#!/usr/bin/env bash


#==============================================================================#
# Script Name: 03b_read_alignment.sh
#
# Last updated: 14/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   Align reads to the Homo sapiens reference genome.
#   []
#
# Usage:
#   Execute from script directory using:
#     sbatch 03b_read_alignment.sh
#
# Software:
#   bowtie v1.3.1
#
# VERSION: 1.0
#
# Command Descriptions:
#   See README.md located in directory containing script for detailed 
#   command descriptions and other useful information:
#   - Inputs
#   - Outputs
#   - Software
#   - File sizes
#   
#==============================================================================#


#==========================#
# SLURM SUBMISSION      
#==========================#
#   SLURM Submission: Defines the SLRUM parameters required to execute
#                     all steps of this script

#SBATCH --partition=defq                          # Edit for desired cluster: <cluster> = "defq", "shortq" (example names)
#SBATCH --nodes=1                                      # Number of nodes
#SBATCH --ntasks=1                                     # Number of tasks 
#SBATCH --cpus-per-task=16                              # Number of cores
#SBATCH --mem=32G                                      # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=01:00:00                                # Run time limit (hh:mm:ss)
#SBATCH --job-name=03b_read_alignment                             # Name assigned to job allocation
#SBATCH --output=../../logs/03_alignment/03b_read_alignment/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/03_alignment/03b_read_alignment/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

# Email notifications for SLURM events (optional - uncomment and edit if desired)
# #SBATCH --mail-type=<type> # <type> = "BEGIN", "END", "FAIL", "ALL"
# #SBATCH --mail-user=<user> # <user> = user@exmail.nottingham.ac.uk (example email address)

# Exit immediately if:
# - Command finishes with non-zero status
# - Pipeline error
set -eo pipefail



#==========================#
# LOG HANDLING
#==========================#
#   Log Handling: Create unique log file for each script execution
#                 to track script progress and key information

# Define log directory
LOG_DIR=$(realpath "../../logs/03_alignment/03b_read_alignment")

# Verify/create log directory
mkdir -p "${LOG_DIR}"

# Define log file
LOG_ID=1 # Unique identifier
LOG="${LOG_DIR}/${SLURM_JOB_NAME}_$(date '+%y-%m-%d')_log_${LOG_ID}.txt" # Includes script name and date

# Set unique log file name
# Loop to increment log ID identifier by 1 until unique ID found
while [[ -f "$LOG" ]]; do
    LOG_ID=$((LOG_ID + 1)) # Increase identifier by 1
    LOG="${LOG_DIR}/${SLURM_JOB_NAME}_$(date '+%y-%m-%d')_log_${LOG_ID}.txt" # Save log file with unique identifier
done 


# Script initialisation message
echo "<------------------------------------------------------->" >> "$LOG"
echo "Initialising script: $SLURM_JOB_NAME.sh" >> "$LOG"
echo "Date initialised" "$(date)" >> "$LOG"


#==========================#
# SETUP NAVIGATION
#==========================#
#   Setup Navigation: Defines the project root to enable 
#                     downstream relative path navigation

# Initiation message
echo "==> Setting up in-script navigation" >> "$LOG"

#   1. Set project root (2 levels above script directory)
PROJECT_ROOT="$(realpath "${SLURM_SUBMIT_DIR}/../..")" 
echo "Project Root: $PROJECT_ROOT" >> "$LOG"

#   2. Extract script prefix for methods subsection
SCRIPT_PREFIX="${SLURM_JOB_NAME%%[!0-9]*}"
echo "Script Prefix: $SCRIPT_PREFIX" >> "$LOG"

#   3. Define this script
SCRIPT=${PROJECT_ROOT}/scripts/${SCRIPT_PREFIX}*/${SLURM_JOB_NAME}.sh
echo "Script Name: $SLURM_JOB_NAME.sh" >> "$LOG"

#   4. Ensure script is called from directory containing script
#      Exit status 1 if script not found
if [ ! -f $SCRIPT ]; then
   echo "Error: script must be called from the directory containing $SLURM_JOB_NAME.sh" >&2
   echo "       Use \"cd <PROJECT_ROOT>/scripts/<SCRIPT_DIRECTORY>\" then re-run the script" >&2
   exit 1
fi

# Completion message
echo "==> Setting up in-script navigation: Finished" >> "$LOG"

#==========================#
# SETUP ENVIRONMENT      
#==========================#
#   Setup Environment:  Activates the appropriate Conda environment
#                       which contains the relevant software

# Initiation message
echo "==> Setting up environment" >> "$LOG"

# Define Conda environment
CONDA_ENV="L4137_01_Protease_bowtie"

# Activate Conda environment:
#   1) Ensure bash profile exists (exit status 1 if profile not found)
#   2) Source bash profile
#   3) Activate Conda environment
if [ -f "$HOME/.bash_profile" ]; then
    source "$HOME/.bash_profile"
else
    echo "Error: Bash profile not found." >&2
    echo "Follow Conda installation instructions detailed <PROJECT_ROOT>/installation.md" >&2
    exit 1
fi
conda activate "$CONDA_ENV"
echo "Conda Envrionment:" "$CONDA_ENV" >> "$LOG"

# Completion message
echo "==> Setting up environment: Finished" >> "$LOG"

#==========================#
# FUNCTIONS      
#==========================#
#   Functions: Contains all user-defined functions called
#             within this script

# FUNCTION 1
# DESCRIPTION: Print the current time and date as "yyyy-mm-dd hh:mm:ss"
#              as this is cleaner than the default output from "date"            
# ARGS: None
timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}



#==========================#
# INPUTS      
#==========================#
#   Inputs: Contains all user-defined input files 
#           and directories called within this script

###==== Files ====###
# Initiation message
echo "==> Setting input files" >> "$LOG"

# Reference genome FASTA file
REF_GENOME="${PROJECT_ROOT}/data/reference/GRCh38.p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Mature miRNA FASTA file (miRBase)
MAT_MIRNA="${PROJECT_ROOT}/data/reference/miRNA/mature/mature_hsa_excl_whitespace.fa"

# Hairpin miRNA FASTA file (miRBase)
HAIR_MIRNA="${PROJECT_ROOT}/data/reference/miRNA/hairpin/hairpin_hsa_excl_whitespace.fa"

# All FASTQ files containing the trimmed NGS sequencing output
TRIMMED_READS=(${PROJECT_ROOT}/results/methods_sections/02_quality_control/02b_trimming/cutadapt/*fastq.gz)

echo "Input file(s):" >> "$LOG"
echo "$REF_GENOME" >> "$LOG"
echo "${TRIMMED_READS[@]}" >> "$LOG"

# Completion message
echo "==> Setting input files: Finished" >> "$LOG"


###==== Directories ====###

# Initiation message
echo "==> Setting input directories" >> "$LOG"

# No input directories required
echo "No input directories required" >> "$LOG"

# Completion message
echo "==> Setting input directories: Finished" >> "$LOG"



#==========================#
# OUTPUTS      
#==========================#
#   Outputs: Contains all user-defined output files 
#            and directories called within this script

###==== Files ====###
# Initiation message
echo "==> Setting output files" >> "$LOG"

# [] sam, bam, sorted bam

# Completion message
echo "==> Setting output files: Finished" >> "$LOG"



###==== Directories ====###
# Initiation message
echo "==> Setting output directories" >> "$LOG"

# Directories for reference bowtie index files 
#   1. Homo sapiens reference genome
GENOME_INDEX_DIR="${PROJECT_ROOT}/data/reference/GRCh38.p14" 

#   2. Mature miRNA
MAT_MIRNA_INDEX_DIR="${PROJECT_ROOT}/data/reference/miRNA/mature" 

#   3. Hairpin miRNA
HAIR_MIRNA_INDEX_DIR="${PROJECT_ROOT}/data/reference/miRNA/hairpin" 

# Directory for all output quality reports
OUTDIR_BOWTIE="${PROJECT_ROOT}/results/methods_sections/03_alignment/03b_read_alignment/bowtie/GRCh38_p14"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_BOWTIE" "$GENOME_INDEX_DIR" "$MAT_MIRNA_INDEX_DIR" "$HAIR_MIRNA_INDEX_DIR")
echo "Output directories required:" >> "$LOG"
echo "${DIR_LIST[@]}" >> $LOG

# Completion message
echo "==> Setting output directories: Finished" >> "$LOG"


#==========================#
# DIRECTORY CREATION      
#==========================#
#   Directory Creation: Ensures output directories exist;
#                       creates the directory, if not

# Initiation message
echo "==> Verifying/creating output directories" >> "$LOG"
 
# For each required output directory:
#   1) Check if directory already exists 
#   2) Creates directory if it doesn't exist
#   3) Prints completion message defining action taken
for directory in "${DIR_LIST[@]}"; do
    if [ -d "${directory}" ]; then
        echo "Directory exists:" "${directory}" >> "$LOG"
    else
        mkdir -p "${directory}"
        echo "Directory created:" "${directory}" >> "$LOG"
    fi
done

# Completion message
echo "==> Verifying/creating output directories: Finished" >> "$LOG"


#==========================#
# CONFIGURATION      
#==========================#
#   Configuration: Sets key variables for downstream use
#                  (e.g., quality thresholds)

# Initiation message
echo "==> Setting parameters" >> "$LOG"

# No specific configuration necessary
echo "No user-defined configuration necessary" >> "$LOG"

# Completion message
echo "==> Setting parameters: Finished" >> "$LOG"


#==========================#
# MAIN SCRIPT      
#==========================#
#   Main Script: Executes the main body of code which produces 
#                the outputs of the script


# Initiation message
echo "$(timestamp)" "==> Initiating bowtie" >> "$LOG"

# Build bowtie index
#   1. Homo sapiens reference genome
#   Index directory prefix = "GRCh38_p14"
bowtie-build "$REF_GENOME" "$GENOME_INDEX_DIR/GRCh38_genome"

echo "$(timestamp)" "==> Created reference genome index" >> "$LOG"

#   2. Mature miRNA
#   Index directory prefix = "mature_hsa_mirbase"
bowtie-build "$MAT_MIRNA" "$MAT_MIRNA_INDEX_DIR/mature_hsa_mirbase"

echo "$(timestamp)" "==> Created mature miRNA index" >> "$LOG"

#   3. Hairpin miRNA
#   Index directory prefix = "hairpin_hsa_mirbase"
bowtie-build "$HAIR_MIRNA" "$HAIR_MIRNA_INDEX_DIR/hairpin_hsa_mirbase"

echo "$(timestamp)" "==> Created hairpin miRNA index" >> "$LOG"


# # Align reads to Homo sapiens reference genome
# for sample in "${TRIMMED_READS[@]}"; do

#     # Extract sample name
#     SAMPLE_NAME=$(basename "$sample" .fastq.gz)

# # Perform alignment
#     bowtie \
#         -a \
#         --best \
#         --strata \
#         -v 2 \
#         --norc \
#         -x "$GENOME_INDEX_DIR" \
#         -p $SLURM_CPUS_PER_TASK \
#         -q "$sample" \
#         -S "$OUTDIR_BOWTIE/${SAMPLE_NAME}.sam"

# # Alignment statistics prep: convert SAM to BAM, sort and index

#     # Convert file
#     samtools view -bS "$OUTDIR_BOWTIE/${SAMPLE_NAME}.sam" > "$OUTDIR_BOWTIE/${SAMPLE_NAME}.bam"

#     # Sort
#     samtools sort -o "$OUTDIR_BOWTIE/${SAMPLE_NAME}_sorted.bam" \
#     "$OUTDIR_BOWTIE/${SAMPLE_NAME}.bam" 

#     # Index
#     samtools index "$OUTDIR_BOWTIE/${SAMPLE_NAME}_sorted.bam"


# # Calculate alignment statistics

#     # Flagstats
#     samtools flagstat "$OUTDIR_BOWTIE/${SAMPLE_NAME}_sorted.bam" > "$OUTDIR_BOWTIE/${SAMPLE_NAME}_flagstat.txt"

#     # Idxstats
#     samtools idxstats "$OUTDIR_BOWTIE/${SAMPLE_NAME}_sorted.bam" > "$OUTDIR_BOWTIE/${SAMPLE_NAME}_idxstats.txt"

# done


# Completion message
echo "$(timestamp)" "==> Bowtie Finished" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
