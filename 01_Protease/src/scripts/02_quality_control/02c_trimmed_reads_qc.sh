#!/usr/bin/env bash


#==============================================================================#
# Script Name: 02c_trimmed_reads_qc.sh
#
# Last updated: 11/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   - Perform quality control (QC) testing of trimmed FASTQ files.
#   - Examine various quality metrics per sample.
#   - Compile QC results for all samples into one report.
#
# Usage:
#     sbatch <PATH_TO_SCRIPT>/quality_control.sh
#
# Software:
#   FastQC v0.12.1
#   miRTrace v1.0.1
#   MultiQC v1.34
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

#SBATCH --partition=defq                                                     # Edit for desired cluster: <cluster> = "defq", "shortq" (example names)
#SBATCH --nodes=1                                                            # Number of nodes
#SBATCH --ntasks=1                                                           # Number of tasks 
#SBATCH --cpus-per-task=4                                                    # Number of cores
#SBATCH --mem=32G                                                            # Memory allocation ("M" = mb, "G" = gb)
#SBATCH --time=02:00:00                                                      # Run time limit (hh:mm:ss)
#SBATCH --job-name=02c_trimmed_reads_qc                                       # Name assigned to job allocation
#SBATCH --output=../../logs/02_quality_control/02c_trimmed_reads_qc/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/02_quality_control/02c_trimmed_reads_qc/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

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
LOG_DIR=$(realpath "../../logs/02_quality_control/02c_trimmed_reads_qc")

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
CONDA_ENV="L4137_01_Protease"

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
echo "Conda Envrionment: " "$CONDA_ENV" >> "$LOG"

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

# All FASTQ files containing the NGS sequencing output as array
TRIMMED_READS=("${PROJECT_ROOT}"/results/methods_sections/02_quality_control/02b_trimming/cutadapt/*.fastq.gz) 
echo "Input file(s):" >> "$LOG"
echo "${TRIMMED_READS[@]}" >> "$LOG"

# Completion message
echo "==> Setting input files: Finished" >> "$LOG"

###==== Directories ====###

# Initiation message
echo "==> Setting input directories" >> "$LOG"

# Parent directory of all QC reports for individual samples from:
# 1. FastQC
# 2. miRTrace
INDIR_SAMPLE_REPORTS="${PROJECT_ROOT}"/results/methods_sections/02_quality_control/02c_trimmed_reads_qc
echo "Input directory:" >> "$LOG"
echo "${INDIR_SAMPLE_REPORTS}" >> $LOG

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
OUTDIR_FASTQC="${PROJECT_ROOT}/results/methods_sections/02_quality_control/02c_trimmed_reads_qc/fastqc"
OUTDIR_MIRTRACE_TRACE="${PROJECT_ROOT}/results/methods_sections/02_quality_control/02c_trimmed_reads_qc/mirtrace/trace"
OUTDIR_MIRTRACE_QC="${PROJECT_ROOT}/results/methods_sections/02_quality_control/02c_trimmed_reads_qc/mirtrace/qc"
OUTDIR_MULTIQC="${PROJECT_ROOT}/results/methods_sections/02_quality_control/02c_trimmed_reads_qc/multiqc"

# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_FASTQC" "$OUTDIR_MIRTRACE_TRACE" "$OUTDIR_MIRTRACE_QC" "$OUTDIR_MULTIQC")
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
        echo "Directory exists: " "${directory}" >> "$LOG"
    else
        mkdir -p "${directory}"
        echo "Directory created: " "${directory}" >> "$LOG"
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

# Check if FASTQ files exist
# Exit status 1 if files not found
if [ ${#TRIMMED_READS[@]} -eq 0 ]; then
    echo "Error: No FASTQ files found in ${PROJECT_ROOT}/data/raw/"
    echo "Execute 01_data_preparation then re-run this script"
    exit 1
fi


###==== FastQC ====###

# Initiation message
echo "$(timestamp)" "==> Initiating FastQC" >> "$LOG"


# Run FastQC to generate QC reports
fastqc \
"${TRIMMED_READS[@]}" \
-o "$OUTDIR_FASTQC" \
-t $SLURM_CPUS_PER_TASK

# Completion message
echo "$(timestamp)" "==> FastQC Finished" >> "$LOG"



###==== miRTrace ====###

# Initiation message
echo "$(timestamp)" "==> Initiating miRTrace" >> "$LOG"


# Run miRTrace to generate trace QC reports
mirtrace \
trace \
-o "$OUTDIR_MIRTRACE_TRACE" \
-f \
-t $SLURM_CPUS_PER_TASK \
"${TRIMMED_READS[@]}" \

# Run miRTrace to generate QC reports
mirtrace \
qc \
-s hsa \
-o "$OUTDIR_MIRTRACE_QC" \
-f \
-t $SLURM_CPUS_PER_TASK \
"${TRIMMED_READS[@]}" \

# Completion message
echo "$(timestamp)" "==> miRTrace Finished" >> "$LOG"

###==== MultiQC ====###

# Initiation message
echo "$(timestamp)" "==> Initiating MultiQC" >> "$LOG"


# Run MultiQC to merge all individual sample QC reports into one report
multiqc \
-o "$OUTDIR_MULTIQC" \
"$INDIR_SAMPLE_REPORTS"

# Completion message
echo "$(timestamp)" "==> MultiQC Finished" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
