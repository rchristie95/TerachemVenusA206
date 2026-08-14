#!/usr/bin/env bash
# Detach as a user systemd unit. setsid+nohup is NOT sufficient -- the agent
# harness kills its whole process tree on exit, and a setsid child is still a
# descendant. systemd reparents the queue outside that tree, so it survives
# logout, the terminal closing, and the agent exiting.
cd "$(dirname "$0")" || exit 1
systemctl --user reset-failed weekend 2>/dev/null
systemd-run --user --unit=weekend --working-directory="$PWD" ./run_weekend.sh
echo "weekend queue running as user systemd unit 'weekend'"
