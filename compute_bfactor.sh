# !/bin/bash
# Computes lysine B-factors of the PDB structures (dsa.bfactor.generate_bfactors) with parallel workers,
# one tmux window per worker. Workers pick random structures, so they rarely collide.
N_WORKERS=${1:-24}
UPDATE_INTERVAL=${2:-200}
ENV_DIR=/data/shajain/container_setup/miniforge3
LOG_DIR=/home/shajain/DynamicSurfaceAccessibility/src/dsa/data/bfactor/logs
mkdir -p $LOG_DIR
tmux new-session -d -s bfactor 2>/dev/null
for i in $(seq 1 $N_WORKERS); do
    tmux new-window -t bfactor "bash -c 'source $ENV_DIR/etc/profile.d/conda.sh && conda activate dsa && python -u -c \"from dsa.bfactor.generate_bfactors import main; main(update_interval=$UPDATE_INTERVAL)\" 2>&1 | tee $LOG_DIR/worker_$i.log; bash'"
done
