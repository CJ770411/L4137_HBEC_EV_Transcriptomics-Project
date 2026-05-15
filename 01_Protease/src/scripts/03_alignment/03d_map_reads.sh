#!/usr/bin/env bash


#==============================================================================#
# Script Name: 03d_map_reads.sh
#
# Last updated: 14/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   []
#
# Usage:
#   Execute from script directory using:
#     sbatch 03d_map_reads.sh
#
# Software:
#   miRDeep2
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
#SBATCH --mem=20G                                      # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=02:00:00                                # Run time limit (hh:mm:ss)
#SBATCH --job-name=03d_map_reads                             # Name assigned to job allocation
#SBATCH --output=../../logs/03_alignment/03d_map_reads/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/03_alignment/03d_map_reads/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --array=0-14 

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
LOG_DIR=$(realpath "../../logs/03_alignment/03d_map_reads")

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
CONDA_ENV="L4137_01_Protease_mirdeep2"

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

###==== Directories ====###

# Initiation message
echo "==> Setting input files" >> "$LOG"
echo "Input file(s):" >> "$LOG"

# Homo sapiens reference genome bowtie index directory
INDIR_GENOME_INDEX="${PROJECT_ROOT}/data/reference/GRCh38.p14/GRCh38_genome" 
echo "$(timestamp)" "==> Directory containing Homo sapiens genome index files: "$INDIR_GENOME_INDEX >> "$LOG"

# Trimmed reads directory
INDIR_TRIMMED_READS="${PROJECT_ROOT}/results/methods_sections/02_quality_control/02b_trimming/cutadapt"
echo "$(timestamp)" "==> Directory containing trimmed reads: "$INDIR_TRIMMED_READS >> "$LOG"

# Completion message
echo "==> Setting input directories: Finished" >> "$LOG"



###==== Files ====###
# Initiation message
echo "==> Setting input files" >> "$LOG"
echo "Input file(s):" >> "$LOG"



# Load trimmed read files into an array for parallel analysis
mapfile -t TRIMMED_READS < <(find "$INDIR_TRIMMED_READS" -name "*.fastq.gz" | sort)

# Identify sample based on SLURM_ARRAY_TASK_ID
SAMPLE_FILE="${TRIMMED_READS[$SLURM_ARRAY_TASK_ID]}"
echo "$(timestamp)" "==> Sample file: "$SAMPLE_FILE >> "$LOG"

# Extract sample name
SAMPLE_NAME=$(basename "$SAMPLE_FILE" .fastq.gz)
echo "$(timestamp)" "==> Sample name: "$SAMPLE_NAME >> "$LOG"


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

# Outputs from alignment (mapper.pl)
OUTDIR_MAP="${PROJECT_ROOT}/results/methods_sections/03_alignment/03d_map_reads/mapper"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_MAP")
echo "Output directories required:" >> "$LOG"
echo "${DIR_LIST[@]}" >> $LOG

# Completion message
echo "==> Setting output directories: Finished" >> "$LOG"



###==== Files ====###
# Initiation message
echo "==> Setting output files" >> "$LOG"

COLLAPSED_READS="$OUTDIR_MAP/${SAMPLE_NAME}_collapsed.fasta"
MAPPED_READS="$OUTDIR_MAP/${SAMPLE_NAME}_vs_genome_GRCh38.arf" 



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

# Mapping filters:
MIN_LENGTH=18

echo "Mapping Filters:" >> "$LOG" 
echo "MIN_LENGTH=$MIN_LENGTH" >> "$LOG"  

# Completion message
echo "==> Setting parameters: Finished" >> "$LOG"

#==========================#
# MAIN SCRIPT      
#==========================#
#   Main Script: Executes the main body of code which produces 
#                the outputs of the script


# Initiation message
echo "$(timestamp)" "==> Initiating miRDeep2" >> "$LOG"

# Define temporary unzipped FASTQ file for mapper.pl (miRDeep2) compatability
TEMP_TRIMMED_FASTQ=$(basename "$TRIMMED_READS" .gz)

# Create temporary unzipped FASTQ file
gunzip -c "$SAMPLE_FILE" > "$TEMP_TRIMMED_FASTQ"

# Execute mapping
mapper.pl \
"$TEMP_TRIMMED_FASTQ" \
-g hsa \
-l "$MIN_LENGTH" \
-n -h -e -i -j -m \
-k TGGAATTCTCGGGTGCCAAGG \
-s "$COLLAPSED_READS" \
-p "$INDIR_GENOME_INDEX" \
-t "$MAPPED_READS" 

# Remove temporary FASTQ file
rm "$TEMP_TRIMMED_FASTQ"



# Completion message
echo "$(timestamp)" "==> miRDeep2 Finished" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
