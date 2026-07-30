#!/usr/bin/env bash


#==============================================================================#
# Script Name: 05_read_quantification.sh
#
# Last updated: 18/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   - Map reads to known mature miRNA sequences to isolate miRNAs
#   - Quantify miRNA expression to enable differential expression analysis
#
# Usage:
#   Execute from script directory using:
#     sbatch 05_read_quantification.sh
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
#SBATCH --mem=16G                                      # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=01:00:00                                # Run time limit (hh:mm:ss)
#SBATCH --job-name=05_read_quantification                             # Name assigned to job allocation
#SBATCH --output=../../../logs/calu3/05_read_quantification/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../../logs/calu3/05_read_quantification/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --array=0-14                                   # One job per sample (15 samples)

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
LOG_DIR=$(realpath "../../../logs/calu3/05_read_quantification")

# Verify/create log directory
mkdir -p "${LOG_DIR}"

# Define unique log file
#    File name: contains script execution-specific information and date executed
LOG="${LOG_DIR}/$(date '+%y-%m-%d')_slurm-${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" 

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

#   1. Set project root (3 levels above script directory)
PROJECT_ROOT="$(realpath "${SLURM_SUBMIT_DIR}/../../..")" 
echo "Project Root: $PROJECT_ROOT" >> "$LOG"

#   3. Define this script
SCRIPT=${PROJECT_ROOT}/scripts/workflow/calu3/${SLURM_JOB_NAME}.sh
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
echo "Conda environment:" "$CONDA_ENV" >> "$LOG"

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


###==== Directories ====###

# Initiation message
echo "==> Setting input directories" >> "$LOG"

# No input directories required
echo "No input directories required" >> "$LOG"


# Completion message
echo "==> Setting input directories: Finished" >> "$LOG"




###==== Files ====###
# Initiation message
echo "==> Setting input files" >> "$LOG"

# Collapsed read files (outoput from mapper.pl) into an array for parallel analysis
COLLAPSED_FILES="${PROJECT_ROOT}/results/calu3/04_collapse_reads/mapper"

# Load collapsed reads into an array for parallel analysis
mapfile -t COLLAPSED_FILES_ARRAY < <(find "$COLLAPSED_FILES" -name "*_collapsed.fasta" | sort)

# Identify sample based on SLURM_ARRAY_TASK_ID
INFILE_COLLAPSED_READS="${COLLAPSED_FILES_ARRAY[$SLURM_ARRAY_TASK_ID]}"
echo  "==> Sample file: "$INFILE_COLLAPSED_READS >> "$LOG"

# Extract sample name
SAMPLE_NAME=$(basename "$INFILE_COLLAPSED_READS" _collapsed.fasta)
echo  "==> Sample name: "$SAMPLE_NAME >> "$LOG"

# Mature miRNA FASTA file (miRBase)
INFILE_MAT_MIRNA="${PROJECT_ROOT}/data/reference/miRNA/mature/mature_hsa_excl_whitespace.fa"
echo  "==> Mature miRNA: "$INFILE_MAT_MIRNA >> "$LOG"

# Hairpin miRNA FASTA file (miRBase)
INFILE_HAIR_MIRNA="${PROJECT_ROOT}/data/reference/miRNA/hairpin/hairpin_hsa_excl_whitespace.fa"
echo  "==> Hairpin miRNA: "$INFILE_HAIR_MIRNA >> "$LOG"



# Completion message
echo "==> Setting input files: Finished" >> "$LOG"




#==========================#
# OUTPUTS      
#==========================#
#   Outputs: Contains all user-defined output files 
#            and directories called within this script


###==== Directories ====###
# Initiation message
echo "==> Setting output directories" >> "$LOG"

# Directory for quantifier.pl outputs
#       Distinct directory for each sample analysed
OUTDIR_SAMPLE="${PROJECT_ROOT}/results/calu3/05_read_quantification/${SAMPLE_NAME}"

#       Directory containing only count matrices for all samples
OUTDIR_ALL_SAMPLES="${PROJECT_ROOT}/results/calu3/05_read_quantification/all_samples"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_SAMPLE" "$OUTDIR_ALL_SAMPLES")
echo "Output directories required:" >> "$LOG"
echo "${DIR_LIST[@]}" >> $LOG

# Completion message
echo "==> Setting output directories: Finished" >> "$LOG"



###==== Files ====###
# Initiation message
echo "==> Setting output files" >> "$LOG"

# No user-defined output files required
echo "No user-defined output files required" >> "$LOG"

# Completion message
echo "==> Setting output files: Finished" >> "$LOG"


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
echo "$(timestamp)" "==> Initiating miRDeep2 for" "$SAMPLE_NAME" >> "$LOG"


# Navigate to the unique sample output directory
cd "$OUTDIR_SAMPLE"

# Quantify read counts to produce count matrix
quantifier.pl \
-T $SLURM_CPUS_PER_TASK \
-p "$INFILE_HAIR_MIRNA" \
-m "$INFILE_MAT_MIRNA" \
-r "$INFILE_COLLAPSED_READS" \
-t hsa



# Identify mature miRNA count file
#       Note: "all_samples" is default naming convention from miRDeep2.
#              The file only contains data from one sample.
COUNT_MATRIX=$(find ./ -name "miRNAs_expressed_all_samples_*.csv")

# Copy count matrix to directory containing all samples and re-name as sample name
cp $COUNT_MATRIX "${OUTDIR_ALL_SAMPLES}"/${SAMPLE_NAME}_count_matrix.csv


# Completion message
echo "$(timestamp)" "==> miRDeep2 Finished for" "$SAMPLE_NAME" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
