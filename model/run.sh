#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

INPUT_FILE="$1"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Define number of parallel jobs
NUM_JOBS=16

# Create log directory
LOG_DIR="$HOME/parallel_logs"
mkdir -p "$LOG_DIR"

chmod 777 "$LOG_DIR"  # Give full read/write/execute permissions

# Log file for parallel execution details
JOB_LOG="$LOG_DIR/parallel.log"

# Print Start Time
echo "Starting Parallel Execution at $(date)" | tee -a "$JOB_LOG"

# Run commands in parallel using input file
parallel --eta -j "$NUM_JOBS" --joblog "$JOB_LOG" \
    '{}; echo "Job {#} completed at $(date)"' < "$INPUT_FILE"

# Print End Time
echo "Parallel Execution Completed at $(date)" | tee -a "$JOB_LOG"

# Display job summary
echo "Job Summary:"
cat "$JOB_LOG"

