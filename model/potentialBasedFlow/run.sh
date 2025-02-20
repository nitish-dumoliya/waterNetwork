#!/bin/bash

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

# Run commands in parallel
parallel --eta -j "$NUM_JOBS" --joblog "$JOB_LOG" \
    '{}; echo "Job {#} completed at $(date)"' < scriptHeuristic.par

# Print End Time
echo "Parallel Execution Completed at $(date)" | tee -a "$JOB_LOG"

# Display job summary
echo "Job Summary:"
cat "$JOB_LOG"

