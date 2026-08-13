#!/usr/bin/env bash
# Detach the ladder from this session entirely: new session id (setsid), no
# controlling terminal, SIGHUP ignored (nohup). It keeps running when the SSH
# connection drops, the network goes down, or the terminal is closed.
cd "$(dirname "$0")" || exit 1
setsid nohup ./run_ladder.sh > ladder_console.log 2>&1 < /dev/null &
echo "ladder detached, pid $!"
