#!/usr/bin/env bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20g
#SBATCH --time=48:00:00
#SBATCH --job-name=run
#SBATCH --output=../../logs/run/slurm-%x-%j.out


##########################

## Author: Chris Janschke
## Date: 29.07.2026
## Description: Script to execute the study workflow.
## Usage: 
#    Execute from script directory

## Software used:
#  Software defined in constituent scripts.


##########################


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
LOG_DIR=$(realpath "../../logs/run")

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
SCRIPT=${PROJECT_ROOT}/scripts/workflow/${SLURM_JOB_NAME}.sh
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

# Source bash profile to enable conda
source $HOME/.bash_profile





#==========================#
# EXECUTE WORKFLOW - Primary HBECs
#==========================#

# Navigate to script directory
cd "${PROJECT_ROOT}/scripts/workflow/primaryHBEC"


## 1. Raw read trimming

# Script: 01_raw_read_trimming.sh


# Execute script
JOB_01=$(sbatch --parsable 01_raw_read_trimming.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_01 --wrap \
"echo \"\$(date): Completed script: 01_raw_read_trimming.sh\" >> '$LOG'"



## 2. Sequencing Reads Quality Control


# Script: 02_reads_qc.sh

# Execute script
JOB_02=$(sbatch --parsable --dependency=afterok:$JOB_01 02_reads_qc.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_02 --wrap \
"echo \"\$(date): Completed script: 02_reads_qc.sh\" >> '$LOG'"



## 3. Download miRNA reference databases


# Script: 03_download_mirna_references.sh

# Execute script
JOB_03=$(sbatch --parsable --dependency=afterok:$JOB_02 03_download_mirna_references.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_03 --wrap \
"echo \"\$(date): Completed script: 03_download_mirna_references.sh\" >> '$LOG'"



## 4. Collapse sequencing reads


# Script: 04_collapse_reads.sh

# Execute script
JOB_04=$(sbatch --parsable --dependency=afterok:$JOB_03 04_collapse_reads.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_04 --wrap \
"echo \"\$(date): Completed script: 04_collapse_reads.sh\" >> '$LOG'"



## 5. Quantify miRNA expression


# Script: 05_read_quantification.sh

# Execute script
JOB_05=$(sbatch --parsable --dependency=afterok:$JOB_04 05_read_quantification.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_05 --wrap \
"echo \"\$(date): Completed script: 05_read_quantification.sh\" >> '$LOG'"




#==========================#
# EXECUTE WORKFLOW - Calu-3
#==========================#

# Navigate to script directory
cd "${PROJECT_ROOT}/scripts/workflow/calu3"


## 6. Raw read trimming

# Script: 01_raw_read_trimming.sh


# Execute script
JOB_06=$(sbatch --parsable --dependency=afterok:$JOB_05 01_raw_read_trimming.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_06 --wrap \
"echo \"\$(date): Completed script: 01_raw_read_trimming.sh\" >> '$LOG'"



## 7. Sequencing Reads Quality Control


# Script: 02_reads_qc.sh

# Execute script
JOB_07=$(sbatch --parsable --dependency=afterok:$JOB_06 02_reads_qc.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_07 --wrap \
"echo \"\$(date): Completed script: 02_reads_qc.sh\" >> '$LOG'"



## 8. Download miRNA reference databases


# Script: 03_download_mirna_references.sh

# Execute script
JOB_08=$(sbatch --parsable --dependency=afterok:$JOB_07 03_download_mirna_references.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_08 --wrap \
"echo \"\$(date): Completed script: 03_download_mirna_references.sh\" >> '$LOG'"



## 9. Collapse sequencing reads


# Script: 04_collapse_reads.sh

# Execute script
JOB_09=$(sbatch --parsable --dependency=afterok:$JOB_08 04_collapse_reads.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_09 --wrap \
"echo \"\$(date): Completed script: 04_collapse_reads.sh\" >> '$LOG'"



## 10. Quantify miRNA expression


# Script: 05_read_quantification.sh

# Execute script
JOB_10=$(sbatch --parsable --dependency=afterok:$JOB_09 05_read_quantification.sh)

# Record completion of script in log file
sbatch --dependency=afterok:$JOB_10 --wrap \
"echo \"\$(date): Completed script: 05_read_quantification.sh\" >> '$LOG'"


