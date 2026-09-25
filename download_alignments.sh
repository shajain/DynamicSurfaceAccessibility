# !/bin/bash
# source /home/shajain/miniconda3/etc/profile.d/conda.sh   # add this line
# conda activate dsa
for i in {1..5}; do
    tmux new-window -t download_alignments "bash -c 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate dsa && python -m dsa.binding_interfaces.structure_uniprot_alignment.pdbe_downloader; bash'"
done
