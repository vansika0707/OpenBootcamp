#!/bin/bash

# Activate conda environment (edit with your env name)
# source activate your_env_name

# Run analysis scripts
python analyze_eukaryotes.py > euk_stats.txt
python find_orfs.py > orf_results.txt
