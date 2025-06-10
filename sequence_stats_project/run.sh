conda create -n seq_env python=3.10 -y
conda activate seq_env
conda install -c conda-forge biopython -y
python analyze_sequences.py > results.txt
