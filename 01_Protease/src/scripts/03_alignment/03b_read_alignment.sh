#!/usr/bin/env bash


#==============================================================================#
# Script Name: 03b_read_alignment.sh
#
# Last updated: 13/05/2026 (dd/mm/yyyy)
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
#SBATCH --cpus-per-task=4                              # Number of cores
#SBATCH --mem=16G                                      # Memory allocation ("M" = mb, "G" = gb)
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
LOG_DIR=../../logs/03_alignment/03b_read_alignment

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
echo "Conda Environemnt:" "$CONDA_ENV" >> "$LOG"

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

# Directory for Homo sapiens reference genome index files
INDEX_DIR="${PROJECT_ROOT}/data/reference/GRCh38.p14/index"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$INDEX_DIR")
echo "Output directories required:" >> "$LOG"
echo "${DIR_LIST[@]}" >> $LOG

# Completion message
echo "==> Setting output files: Finished" >> "$LOG"



###==== Directories ====###
# Initiation message
echo "==> Setting output directories" >> "$LOG"

# Directory for all output quality reports
OUTDIR_BOWTIE="${PROJECT_ROOT}/results/methods_sections/03_alignment/03b_read_alignment/bowtie/GRCh38_p14"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_BOWTIE")
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
#   Index directory prefix = "GRCh38_p14"
bowtie-build "$REF_GENOME" "$INDEX_DIR"

# Align reads to Homo sapiens reference genome
for sample in "${TRIMMED_READS[@]}"; do

    SAMPLE_NAME=$(basename "$sample" .fastq.gz)

    bowtie \
        -a \
        --best \
        --strata \
        -v 2 \
        --nofw \
        -x "$INDEX_DIR" \
        -p $SLURM_CPUS_PER_TASK \
        -q "$sample" \
        -S "$OUTDIR_BOWTIE/${SAMPLE_NAME}.sam"

done

# [] Example commands
bowtie -v 1 -m 1 --best --strata --norc -l 20 \
-x ~/Genome_Index/mirbase_bowtie_index/mirbase_hsa \
-p \
-q ~/miRNA/trimmed/SRR1759248/SRR1759248_trimmed.fq.gz \
-S ~/miRNA/aligned/SRR1759248/SRR1759248_trimmed.sam

bowtie -v 0 -a --best --strata reference_index trimmed_reads.fastq > aligned_reads.sam

bowtie -v 2 -k 1 --best --strata \
  grch38_index \
  sample.fastq \
  sample_genome_mapped.bwt



# Completion message
echo "$(timestamp)" "==> Bowtie Finished" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
