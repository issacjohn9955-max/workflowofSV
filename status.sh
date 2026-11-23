#!/usr/bin/env bash
set -euo pipefail
out=$1
python3 - "$out" <<'PY'
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
status=root/'progress.json'
if status.exists():
    state=json.loads(status.read_text())
    n=len(state['completed'])
    print(f"{state['status']} | {n}/{state['total']} ({100*n/state['total']:.0f}%) | {state['stage'] or 'all stages'}")
    print("Completed: "+", ".join(state['completed']))
else:
    print("Not started")
PY
if [[ -f "$out/worker.log" ]]; then tail -8 "$out/worker.log"; fi
