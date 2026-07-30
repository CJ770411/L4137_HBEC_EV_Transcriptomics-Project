#!/usr/bin/env bash


#==============================================================================#
# Script Name: 03_download_mirna_references.sh
#
# Last updated: 13/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   - Download reference files necessary to create bowtie index and perform alignment. 
#     References to download:
#       2. (miRBase) Mature miRNA sequences
#       3. (miRBase) Hairpin miRNA sequences
#   - Clean up raw reference files for compatability with downstream software.
#   - Extract only Homo sapiens data from miRBase files.
#
# Usage:
#   Execute from script directory using:
#     sbatch 03_download_mirna_references.sh
#
# Software:
#   miRDeep2 v2.0.1.3
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
#SBATCH --cpus-per-task=1                              # Number of cores
#SBATCH --mem=8G                                       # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=01:00:00                                # Run time limit (hh:mm:ss)
#SBATCH --job-name=03_download_mirna_references                             # Name assigned to job allocation
#SBATCH --output=../../../logs/primaryHBEC/03_download_mirna_references/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../../logs/primaryHBEC/03_download_mirna_references/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

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
LOG_DIR=$(realpath "../../../logs/primaryHBEC/03_download_mirna_references")

# Verify/create log directory
mkdir -p "${LOG_DIR}"

# Define unique log file
#    File name: contains script execution-specific information and date executed
LOG="${LOG_DIR}/$(date '+%y-%m-%d')_slurm-${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log" 


# Script initialisation message
echo "<------------------------------------------------------->" >> "$LOG"
echo "Initialising script: $SLURM_JOB_NAME.sh" >> "$LOG"
echo "Date initialised:" "$(date)" >> "$LOG"


#==========================#
# SETUP NAVIGATION
#==========================#
#   Setup Navigation: Defines the project root to enable 
#                     downstream relative path navigation

# Initiation message
echo "==> Setting up in-script navigation" >> "$LOG"

#   1. Set project root (3 levels above script directory)
PROJECT_ROOT="$(realpath "${SLURM_SUBMIT_DIR}/../../..")" 
echo "Project Root: $PROJECT_ROOT" >> "$LOG"

#   3. Define this script
SCRIPT=${PROJECT_ROOT}/scripts/workflow/primaryHBEC/${SLURM_JOB_NAME}.sh
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
CONDA_ENV="L4137_mirdeep2"

# Activate Conda environment:
#   1) Ensure bash profile exists (exit status 1 if profile not found)
#   2) Source bash profile
#   3) Activate Conda environment
if [ -f "$HOME/.bash_profile" ]; then
    source "$HOME/.bash_profile"
else
    echo "Error: Bash profile not found." >&2
    exit 1
fi
conda activate "$CONDA_ENV"
echo "Conda Environment: " "$CONDA_ENV" >> "$LOG"

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

# No input files required
echo "No input files required" >> "$LOG"

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

# Directories for reference data files
OUTDIR_MIRNA="${PROJECT_ROOT}/data/reference/miRNA"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_MIRNA/mature" "$OUTDIR_MIRNA/hairpin")
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


# Define the types of miRNA to be downloaded from miRBase
MIRNA_LIST=("mature" "hairpin")

# Download and clean miRNA reference files for 'mature' and 'hairpin'
for mirna_form in "${MIRNA_LIST[@]}"; do

    # Ensure output directory exists; if not, creates it
    mkdir -p "${OUTDIR_MIRNA}/${mirna_form}"

    # Navigate to reference directory
    cd "${OUTDIR_MIRNA}/${mirna_form}"

    # Initiation message
    echo "$(timestamp)" "==> Downloading $mirna_form miRNA sequences (all species) from miRBase" >> "$LOG"

    # Download FASTA reference file (not zipped)
    wget https://www.mirbase.org/download/${mirna_form}.fa

    # Extract Homo sapiens (hsa) miRNA sequences
    awk '/^>/ {keep = ($0 ~ /hsa/)} keep' ${mirna_form}.fa > ${mirna_form}_hsa_inc_whitespace.fa

    # Remove white-space from FASTA file for downstream compatbility with miRDeep2 (alignment software)
    remove_white_space_in_id.pl ${mirna_form}_hsa_inc_whitespace.fa > ${mirna_form}_hsa_excl_whitespace.fa

    # Completion message
    echo "$(timestamp)" "==> Finished downloading $mirna_form miRNA sequences (all species) from miRBase" >> "$LOG"

done




# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed:" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
