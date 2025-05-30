#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

INPUT_FILE="$1"
if [ ! -f "$INPUT_FILE" ]; then
    echo "File '$INPUT_FILE' not found"
    exit 1
fi

LOG_DIR="$HOME/parallel_logs"
mkdir -p "$LOG_DIR"
chmod 777 "$LOG_DIR"
JOB_LOG="$LOG_DIR/parallel.log"

echo "Starting Parallel Execution at $(date)" | tee -a "$JOB_LOG"

parallel --eta -j 14 --joblog "$JOB_LOG" '{}; echo Job {#} done at $(date)' < "$INPUT_FILE"


echo "Parallel Execution Completed at $(date)" | tee -a "$JOB_LOG"

# Display job summary
echo "Job Summary:"
cat "$JOB_LOG"
