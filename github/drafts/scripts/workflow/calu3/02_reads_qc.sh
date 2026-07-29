#!/usr/bin/env bash


#==============================================================================#
# Script Name: 02_reads_qc.sh
#
# Last updated: 11/05/2026 (dd/mm/yyyy)
#
# Purpose:
#   - Perform quality control (QC) testing of trimmed FASTQ files.
#   - Examine various quality metrics per sample.
#   - Compile QC results for all samples into one report.
#
# Usage:
# Usage:
#   Execute from script directory using:
#     sbatch 02_reads_qc.sh
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
#SBATCH --job-name=02_reads_qc                                       # Name assigned to job allocation
#SBATCH --output=../../logs/analysis_pipeline/run_03/02_reads_qc/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/analysis_pipeline/run_03/02_reads_qc/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

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
LOG_DIR=$(realpath "../../logs/analysis_pipeline/run_03/02_reads_qc")

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

#   2. Define this script
SCRIPT=${PROJECT_ROOT}/scripts/analysis_pipeline/${SLURM_JOB_NAME}.sh
echo "Script Name: $SLURM_JOB_NAME.sh" >> "$LOG"

#   3. Ensure script is called from directory containing script
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
RAW_READS=(${PROJECT_ROOT}/data/raw/*fastq.gz)
TRIMMED_READS=("${PROJECT_ROOT}"/results/analysis_pipeline/run_03/01_raw_read_trimming/cutadapt/*.fastq.gz) 
echo "Input file(s):" >> "$LOG"

echo "Raw read file(s):" >> "$LOG"
echo "${RAW_READS[@]}" >> "$LOG"

echo "Trimmed read file(s):" >> "$LOG"
echo "${TRIMMED_READS[@]}" >> "$LOG"

# Completion message
echo "==> Setting input files: Finished" >> "$LOG"

###==== Directories ====###

# Initiation message
echo "==> Setting input directories" >> "$LOG"

# Parent directory of all QC reports for individual samples from:
# 1. FastQC
# 2. miRTrace
INDIR_SAMPLE_REPORTS_RAW="${PROJECT_ROOT}"/results/analysis_pipeline/run_03/02_reads_qc/raw
INDIR_SAMPLE_REPORTS_TRIMMED="${PROJECT_ROOT}"/results/analysis_pipeline/run_03/02_reads_qc/trimmed
echo "Input directory:" >> "$LOG"
echo "${INDIR_SAMPLE_REPORTS_RAW}" >> $LOG
echo "${INDIR_SAMPLE_REPORTS_TRIMMED}" >> $LOG

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
OUTDIR_RAW="${PROJECT_ROOT}/results/analysis_pipeline/run_03/02_reads_qc/raw"
OUTDIR_TRIMMED="${PROJECT_ROOT}/results/analysis_pipeline/run_03/02_reads_qc/trimmed"


# Combine directories into array for simultaenous creation later
DIR_LIST=("$OUTDIR_RAW" "$OUTDIR_TRIMMED")
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



# RAW reads
# Create raw directory
mkdir -p "$OUTDIR_RAW"/fastqc

# Run FastQC to generate QC reports
fastqc \
"${RAW_READS[@]}" \
-o "$OUTDIR_RAW"/fastqc \
-t $SLURM_CPUS_PER_TASK


# TRIMMED reads
# Create trimmed directory
mkdir -p "$OUTDIR_TRIMMED"/fastqc

# Run FastQC to generate QC reports
fastqc \
"${TRIMMED_READS[@]}" \
-o "$OUTDIR_TRIMMED"/fastqc \
-t $SLURM_CPUS_PER_TASK


# Completion message
echo "$(timestamp)" "==> FastQC Finished" >> "$LOG"



###==== miRTrace ====###

# Initiation message
echo "$(timestamp)" "==> Initiating miRTrace" >> "$LOG"

# RAW reads
# Create raw directories
mkdir -p "$OUTDIR_RAW"/mirtrace/trace
mkdir -p "$OUTDIR_RAW"/mirtrace/qc

# Run miRTrace to generate trace QC reports
mirtrace \
trace \
-o "$OUTDIR_RAW"/mirtrace/trace \
-f \
-t $SLURM_CPUS_PER_TASK \
"${RAW_READS[@]}" 

# Run miRTrace to generate QC reports
mirtrace \
qc \
-s hsa \
-o "$OUTDIR_RAW"/mirtrace/qc \
-f \
-t $SLURM_CPUS_PER_TASK \
"${RAW_READS[@]}" 


# TRIMMED reads
# Create raw directories
mkdir -p "$OUTDIR_TRIMMED"/mirtrace/trace
mkdir -p "$OUTDIR_TRIMMED"/mirtrace/qc


# Run miRTrace to generate trace QC reports
mirtrace \
trace \
-o "$OUTDIR_TRIMMED"/mirtrace/trace \
-f \
-t $SLURM_CPUS_PER_TASK \
"${TRIMMED_READS[@]}" 

# Run miRTrace to generate QC reports
mirtrace \
qc \
-s hsa \
-o "$OUTDIR_TRIMMED"/mirtrace/qc \
-f \
-t $SLURM_CPUS_PER_TASK \
"${TRIMMED_READS[@]}" 



# Completion message
echo "$(timestamp)" "==> miRTrace Finished" >> "$LOG"



###==== MultiQC ====###

# Initiation message
echo "$(timestamp)" "==> Initiating MultiQC" >> "$LOG"

# RAW reads
# Create raw directory
mkdir -p "$OUTDIR_RAW"/multiqc

# Run MultiQC to merge all individual sample QC reports into one report
multiqc \
-o "$OUTDIR_RAW"/multiqc \
"$INDIR_SAMPLE_REPORTS_RAW"


# TRIMMED reads
# Create raw directory
mkdir -p "$OUTDIR_TRIMMED"/multiqc

# Run MultiQC to merge all individual sample QC reports into one report
multiqc \
-o "$OUTDIR_TRIMMED"/multiqc \
"$INDIR_SAMPLE_REPORTS_TRIMMED"



# Completion message
echo "$(timestamp)" "==> MultiQC Finished" >> "$LOG"



# Script completion message
echo "Completed script: $SLURM_JOB_NAME.sh"  >> "$LOG"
echo "Date completed" "$(date)" >> "$LOG"
echo "<------------------------------------------------------->" >> "$LOG"

#==============================================================================#
# End of script
#==============================================================================#
