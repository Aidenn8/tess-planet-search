#!/bin/bash
# Wait for the download and the current search pass to finish, then run another
# pass that picks up every star downloaded since (finished stars are skipped).
cd "$(dirname "$0")/.."
while pgrep -f "scripts/02_download.py" > /dev/null; do sleep 30; done
while pgrep -f "scripts/04_search.py" > /dev/null; do sleep 30; done
echo "$(date) starting second search pass" >> data/search_main.log
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python -u scripts/04_search.py --workers 4 >> data/search_main.log 2>&1
echo "$(date) second search pass finished" >> data/search_main.log
