#!/bin/bash -l

#$ -l h_rt=12:00:00             # Max runtime
#$ -P bf528         # REQUIRED: Replace with your SCC project
#$ -N bootstrap            # Job name
#$ -j y                         # Merge stdout/stderr
#$ -o bootstrap_log.txt                   # Output log 
#$ -m beas           # Send email at (b)egin and (e)nd
#$ -M yuansh24@bu.edu  # Replace with your actual BU email
#$ -cwd                         # Run job from current working directory
#$ -pe omp 12                    # Request X cores
#$ -l mem_per_core=8G           # memory per core

# --- Setup ---
module load R/4.4.0

# Optional: Print basic info
echo "Job ID: $JOB_ID"
echo "Running on host: $(hostname)"
echo "Using $NSLOTS cores and $((NSLOTS * 8))G memory"
echo "Current working directory: $(pwd)"
echo "Current date and time: $(date)"
echo "R version: $(R --version)"

# --- Run R script ---
Rscript cv.R