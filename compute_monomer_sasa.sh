# !/bin/bash
# source /home/shajain/miniconda3/etc/profile.d/conda.sh   # add this line
# conda activate dsa
for i in {1..20}; do
    tmux new-window -t monomer_sasa "bash -c 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate dsa && python -m dsa.sasa.compute_sasa; bash'"
done
