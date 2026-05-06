#!/usr/bin/env bash


#==============================================================================#
# Script Name: [].sh
#
# Last updated: []/[]/[] (dd/mm/yyyy)
#
# Purpose:
#   []
#
# Usage:
#   Execute from script directory using:
#     sbatch [].sh
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
#SBATCH --partition=<cluster>                          # Edit for desired cluster: <cluster> = "defq", "shortq" (example names)
#SBATCH --nodes=1                                      # Number of nodes
#SBATCH --ntasks=1                                     # Number of tasks 
#SBATCH --cpus-per-task=1                              # Number of cores
#SBATCH --mem=16g                                      # Memory allocation ("m" = mb, "g" = gb)
#SBATCH --time=01:00:00                                # Run time limit (hh:mm:ss)
#SBATCH --job-name=$0                                  # Name assigned to job allocation
#SBATCH --output=../../logs/$0/slurm-%x-%j.out               # Standard output log file ("%x" is replaced with job name, "%j" is replaced with job ID)
#SBATCH --error=../../logs/$0/slurm-%x-%j.err                # Standard error log file ("%x" is replaced with job name, "%j" is replaced with job ID)

# Email notifications for SLURM events (optional - uncomment and edit if desired)
# #SBATCH --mail-type=<type> # <type> = "BEGIN", "END", "FAIL", "ALL"
# #SBATCH --mail-user=<user> # <user> = user@exmail.nottingham.ac.uk (example email address)

# Exit immediately if:
# - Command finishes with non-zero status
# - Unset variable
# - Pipeline error
set -euo pipefail

#==========================#
# SETUP ENVIRONMENT      
#==========================#

# Define Conda environment
CONDA_ENV="[]"

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
timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


#==========================#
# CONFIGURATION      
#==========================#



### Inputs ###


### Outputs ###




#==========================#
# DIRECTORY CREATION      
#==========================#
# For each required output directory:
#   1) Check if directory already exists 
#   2) Creates directory if it doesn't exist
#   3) Prints completion message defining action taken

DIR1='sdfsdf'
DIR2='sdfsdf'
DIR3='sdfsdf'
DIR4='sdfsdf'

DIR_LIST="$DIR1" "$DIR2" "$DIR3" "$DIR4"

for directory in "$DIR_LIST"; do
    if [ -d "${directory}" ]; then
    printf "\n%s: Directory exists: %s\n\n" "$(timestamp)" "${directory}"
else
    mkdir -p "${directory}"
    printf "\n%s: Directory created: %s\n\n" "$(timestamp)" "${directory}"
    fi
done

#[check if this loop works]


# Directory: []
if [ -d "$DIR" ]; then
    printf "\n%s: Directory exists: %s\n\n" "$(timestamp)" "$DIR"
else
    mkdir -p "$DIR"
    printf "\n%s: Directory created: %s\n\n" "$(timestamp)" "$DIR"
fi

#[create directory dictionary and create everything with a for loop]




#==========================#
# MAIN SCRIPT      
#==========================#



#==============================================================================#
# End of script
#==============================================================================#
