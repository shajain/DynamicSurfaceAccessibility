# !/bin/bash
# Computes SASA of the PDB biological assemblies (dsa.sasa.generate_sasa_multimers) with parallel workers,
# one tmux window per worker. Workers pick random structures, so they rarely collide.
N_WORKERS=${1:-24}
UPDATE_INTERVAL=${2:-200}
ENV_DIR=/data/shajain/container_setup/miniforge3
LOG_DIR=/home/shajain/DynamicSurfaceAccessibility/src/dsa/data/sasa/multimer_logs
mkdir -p $LOG_DIR
tmux new-session -d -s multimer_sasa 2>/dev/null
for i in $(seq 1 $N_WORKERS); do
    tmux new-window -t multimer_sasa "bash -c 'source $ENV_DIR/etc/profile.d/conda.sh && conda activate dsa && python -u -c \"from dsa.sasa.generate_sasa_multimers import main; main(update_interval=$UPDATE_INTERVAL)\" 2>&1 | tee $LOG_DIR/worker_$i.log; bash'"
done
