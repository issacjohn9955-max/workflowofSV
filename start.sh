#!/usr/bin/env bash
set -euo pipefail
if (($# != 2)); then printf 'Usage: bash start.sh CONFIG.json OUTPUT_DIRECTORY\n' >&2; exit 2; fi
root=$(cd -- "$(dirname -- "$0")" && pwd)
config=$(realpath -e -- "$1")
out=$(realpath -m -- "$2")
mkdir -p -- "$out"
export SV_PYTHON=${SV_PYTHON:-python3}
export OPENBLAS_NUM_THREADS=${SV_THREADS:-1}
export OMP_NUM_THREADS=$OPENBLAS_NUM_THREADS MKL_NUM_THREADS=$OPENBLAS_NUM_THREADS
setsid nohup bash "$root/scripts/worker.sh" "$config" "$out" </dev/null >>"$out/worker.log" 2>&1 &
if command -v tmux >/dev/null && command -v watch >/dev/null; then
    socket="svpc_$(printf '%s' "$out" | sha256sum | cut -c1-12)"
    printf -v monitor 'watch -n 8 -t bash %q %q' "$root/status.sh" "$out"
    tmux -L "$socket" new-session -d -s monitor "$monitor" 2>/dev/null ||
        tmux -L "$socket" has-session -t monitor
    printf 'Monitor: tmux -L %s attach -t monitor\n' "$socket"
fi
printf 'Log: %s/worker.log\n' "$out"
