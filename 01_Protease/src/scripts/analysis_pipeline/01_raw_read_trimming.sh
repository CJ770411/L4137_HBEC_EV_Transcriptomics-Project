#!/usr/bin/env bash


#==============================================================================#
# Script Name: 01_raw_read_trimming.sh
#
# Last updated: 11/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   - Isolate miRNA reads
#   - Remove adapter contamination (PCR by-product)
#
# Usage:
#   Execute from script directory using:
#     sbatch 01_raw_read_trimming.sh
#
# Software:
#   cutadapt v5.2
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

#SBATCH --partition=defq                               # Edit for desired cluster: <cluster> = "defq", "shortq" (example names)
#SBATCH --nodes=1                                      # Number of nodes
#SBATCH --ntasks=1                                     # Number of tasks 
#SBATCH --cpus-per-task=16                             # Number of cores
#SBATCH --mem=32G                                      # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=02:00:00                                # Run time limit (hh:mm:ss)
#SBATCH --job-name=01_raw_read_trimming                             # Name assigned to job allocation
#SBATCH --output=../../logs/analysis_pipeline/run_01/01_raw_read_trimming/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/analysis_pipeline/run_01/01_raw_read_trimming/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

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
LOG_DIR=$(realpath "../../logs/analysis_pipeline/run_01/01_raw_read_trimming")

# Verify/create log directory
mkdir -p "${LOG_DIR}"

# Define unique log file
#    File name: contains script execution-specific information and date executed
LOG="${LOG_DIR}/$(date '+%y-%m-%d')_slurm-${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log" 


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

#   3. Define this script
SCRIPT=${PROJECT_ROOT}/scripts/analysis_pipeline/${SLURM_JOB_NAME}.sh
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
CONDA_ENV="L4137_01_Protease_cutadapt"

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

# All FASTQ files containing the NGS sequencing output
RAW_READS=(${PROJECT_ROOT}/data/raw/*fastq.gz)
echo "Input file(s):" >> "$LOG"
echo "${RAW_READS[@]}" >> "$LOG"

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

# No user-defined output files required
echo "No user-defined output files required" >> "$LOG"

# Completion message
echo "==> Setting output files: Finished" >> "$LOG"


###==== Directories ====###
# Initiation message
echo "==> Setting output directories" >> "$LOG"

# Directory for all output quality reports
OUTDIR_CUTADAPT="${PROJECT_ROOT}/results/analysis_pipeline/run_01/01_raw_read_trimming/cutadapt"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_CUTADAPT")
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


# Trimming filters:
MIN_LENGTH=17
MIN_QUAL=30

echo "Trimming Filters:" >> "$LOG" 
echo "MIN_LENGTH=$MIN_LENGTH" >> "$LOG" 
echo "MIN_QUAL=$MIN_QUAL" >> "$LOG" 


# Completion message
echo "==> Setting parameters: Finished" >> "$LOG"

#==========================#
# MAIN SCRIPT      
#==========================#
#   Main Script: Executes the main body of code which produces 
#                the outputs of the script


# Initiation message
echo "$(timestamp)" "==> Initiating cutadapt" >> "$LOG"

# Run cutadapt to trim raw read data for each sample
for sample in "${RAW_READS[@]}"; do

    SAMPLE_NAME=$(basename "$sample" .fastq.gz)

    cutadapt \
        -j "$SLURM_CPUS_PER_TASK" \
        -a Realseq3P=TGGAATTCTC \
        -u1 \
        --discard-untrimmed \
        -m "$MIN_LENGTH" \
        -q "$MIN_QUAL" \
        -o "${OUTDIR_CUTADAPT}/${SAMPLE_NAME}_trimmed.fastq.gz" \
        "$sample"

done

# Completion message
echo "$(timestamp)" "==> cutadapt Finished" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
