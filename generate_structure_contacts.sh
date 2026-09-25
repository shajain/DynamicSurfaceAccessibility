# !/bin/bash
# source /home/shajain/miniconda3/etc/profile.d/conda.sh   # add this line
# conda activate dsa
for i in {1..5}; do
    #tmux new-window -t generate_afold_contacts "bash -c 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate dsa && python -m dsa.
    #do I need to change anythong now that I am running on container with mamba? If so suggest a command to do so.  
    tmux new-window -t afold_contacts "bash -c 'source /home/shajain/miniforge3/etc/profile.d/conda.sh && conda activate dsa && python -m dsa.binding_interfaces.contacts.generate_contacts; bash'"
done
