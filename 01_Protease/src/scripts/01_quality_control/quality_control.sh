#!/usr/bin/env bash


#==============================================================================#
# Script Name: quality_control.sh
#
# Last updated: 07/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   Perform quality control (QC) testing of raw read FASTQ files.
#
# Usage:
#     sbatch <PATH_TO_SCRIPT>/quality_control.sh
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
#SBATCH --partition=defq                                                     # Edit for desired cluster: <cluster> = "defq", "shortq" (example names)
#SBATCH --nodes=1                                                            # Number of nodes
#SBATCH --ntasks=1                                                           # Number of tasks 
#SBATCH --cpus-per-task=4                                                    # Number of cores
#SBATCH --mem=32G                                                            # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=02:00:00                                                      # Run time limit (hh:mm:ss)
#SBATCH --job-name=quality_control                                           # Name assigned to job allocation
#SBATCH --output=../../logs/01_quality_control/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/01_quality_control/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

# Email notifications for SLURM events (optional - uncomment and edit if desired)
# #SBATCH --mail-type=<type> # <type> = "BEGIN", "END", "FAIL", "ALL"
# #SBATCH --mail-user=<user> # <user> = user@exmail.nottingham.ac.uk (example email address)

# Exit immediately if:
# - Command finishes with non-zero status
# - Unset variable
# - Pipeline error
set -euo pipefail


#==========================#
# SETUP NAVIGATION
#==========================#
# Define project root to enable relative path navigation

#   1. Identify directory containing this script
SCRIPT_DIR="$(dirname "$(realpath "$0")")"

#   2. Set project root (2 levels above script directory)
PROJECT_ROOT="$(realpath "${SCRIPT_DIR}/../..")"

#   3. Confirm file paths
echo "$(date)"
echo "Initiating script: $SLURM_JOB_NAME"
echo "Project Root: $PROJECT_ROOT"
echo "Script Directory: $SCRIPT_DIR"



#==========================#
# SETUP ENVIRONMENT      
#==========================#

# Define Conda environment
CONDA_ENV="L4137_01_Protease"

# Activate Conda environment:
#   1) Ensure bash profile exists (exit status 1 if profile not found)
#   2) Source bash profile
#   3) Activate Conda environment
if [ -f "$HOME/.bash_profile" ]; then
    source "$HOME/.bash_profile"
else
    echo "Error: Bash profile not found."
    echo "Follow Conda installation instructions detailed <PROJECT_ROOT>/installation.md"
    exit 1
fi
conda activate "$CONDA_ENV"


#==========================#
# FUNCTIONS      
#==========================#
# DESCRIPTION: Print the current time and date as "yyyy-mm-dd hh:mm:ss"
#              as this is cleaner than the default output from "date"            
# ARGS: None
timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}



#==========================#
# INPUTS      
#==========================#

###==== Files ====###

# All FASTQ files containing the NGS sequencing output
RAW_READS=${PROJECT_ROOT}/data/raw_data/*fastq.gz

###==== Directories ====###
# No input directories required



#==========================#
# OUTPUTS      
#==========================#

###==== Files ====###

# No output files required



###==== Directories ====###

# Directory for all output quality reports
OUTDIR_FASTQC="${PROJECT_ROOT}/results/reports/fastqc"



# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_FASTQC")



#==========================#
# CONFIGURATION      
#==========================#
# No specific configuration necessary





#==========================#
# DIRECTORY CREATION      
#==========================#
# For each required output directory:
#   1) Check if directory already exists 
#   2) Creates directory if it doesn't exist
#   3) Prints completion message defining action taken

for directory in "${DIR_LIST[@]}"; do
    if [ -d "${directory}" ]; then
        printf "\n%s: Directory exists: %s\n\n" "$(timestamp)" "${directory}"
    else
        mkdir -p "${directory}"
        printf "\n%s - Directory created: %s\n\n" "$(timestamp)" "${directory}"
    fi
done




#==========================#
# MAIN SCRIPT      
#==========================#

# Check if FASTQ files exist
# Exit status 1 if files not found
if [ ${#RAW_READS[@]} -eq 0 ]; then
    echo "Error: No FASTQ files found in ${PROJECT_ROOT}/data/raw_data/"
    echo "Execute 01_preprocessing.sh then re-run this script"
    exit 1
fi

# Run FastQC to generate QC reports
fastqc \
$RAW_READS \
-o "$OUTDIR_FASTQC" \
-t $SLURM_CPUS_PER_TASK




#==============================================================================#
# End of script
#==============================================================================#
