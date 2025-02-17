#!/bin/bash

# Define the number of parallel jobs (adjust based on system capability)
NUM_JOBS=8

# Read commands from a file (e.g., scriptHeuristic.par) and execute them in parallel
parallel --eta -j "$NUM_JOBS" < scriptHeuristic.par
