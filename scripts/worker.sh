#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")" && pwd)
config=$1
out=$2
MAX_TRY=3
RETRY_SECONDS=60
ATTEMPT_TIMEOUT=21600
exec 9>"$out/worker.lock"
if ! flock -n 9; then printf '[%s] Another background worker is active\n' "$(date -Iseconds)"; exit 1; fi
for ((attempt=1; attempt<=MAX_TRY; attempt++)); do
    printf '[%s] START attempt %d/%d\n' "$(date -Iseconds)" "$attempt" "$MAX_TRY"
    if timeout --kill-after=60 "$ATTEMPT_TIMEOUT" "$SV_PYTHON" "$root/workflow.py" "$config" "$out"; then
        printf '[%s] COMPLETE\n' "$(date -Iseconds)"
        exit 0
    fi
    printf '[%s] FAILED attempt %d/%d\n' "$(date -Iseconds)" "$attempt" "$MAX_TRY"
    if ((attempt<MAX_TRY)); then sleep "$RETRY_SECONDS"; fi
done
printf '[%s] FAILED after all attempts\n' "$(date -Iseconds)"
exit 1
