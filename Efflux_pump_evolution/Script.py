"""
Parallel Execution Script for Pseudomonas aeruginosa Evolution Simulations

This script manages the parallel execution of multiple evolutionary simulations
across different parameter sets (rates of regulation) and trial repetitions.

It runs the evolution model for different network configurations in parallel
to explore parameter space and generate statistical results.

Author: Dr. Suvam Roy
Institution: Umeå University
"""

import os
import numpy as np

# Base directory for network parameter files
basic_directory = './nets/'

# First we copy the code just for archive and then we do the loop over parameters...

# Parallel execution loop over network parameter sets and trials
# net: Iterates through different network parameter sets (regulation rates)
# trial: Runs multiple independent trials for each parameter set for statistics

# Loop through network parameter sets 0 to 24 (inclusive)
for net in range(0, 25):
    source = str(net)  # Convert to string for file naming
    
    # For each network parameter set, run 25 independent evolutionary trials
    for trial in range(0, 25):
        target = str(trial)  # Convert to string for file naming
        
        # Construct command to run evolution simulation in background
        # Format: python3 evolution.py [target] [source] &
        # - target: Trial ID (0-24) for file naming and random seed variation
        # - source: Network parameter set ID (0-24) specifying regulation rates
        # - &: Runs process in background allowing parallel execution
        cmd = 'python3 evolution2.py ' + target + ' ' + source + ' &'
        
        # Execute the command
        os.system(cmd)
        
        # Print command for monitoring and debugging
        print(cmd)

# Note: This script will launch 25 networks × 25 trials = 625 parallel processes
# The '&' operator allows all processes to run concurrently in the background
# This significantly reduces total computation time compared to serial execution

# Expected output file structure:
# ./tob/[final_target]genome.npy - evolved genome for each simulation
# ./tob/[final_target]ts.npy - time series data for each simulation
# where final_target = (source * 25 + target) creates unique identifiers
