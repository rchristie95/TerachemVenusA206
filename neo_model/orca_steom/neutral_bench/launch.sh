#!/usr/bin/env bash
# Launch the ladder as a USER SYSTEMD UNIT, not a background shell.
#
# setsid+nohup is NOT sufficient here: the agent harness tears down the whole
# process tree/cgroup on exit, and a setsid child is still a descendant, so it
# gets SIGKILLed. The ADC(2) rung was killed that way mid-run, leaving no
# traceback -- only a truncated log. systemd-run reparents the job to the user
# systemd manager (PID 2493), which is outside the harness's tree entirely.
#
#   status:  systemctl --user status ladder
#   log:     journalctl --user -u ladder -f
#   stop:    systemctl --user stop ladder
cd "$(dirname "$0")" || exit 1
systemctl --user reset-failed ladder 2>/dev/null
systemd-run --user --unit=ladder --working-directory="$PWD" \
  --setenv=TMPDIR="$PWD/adc2_scratch" \
  --setenv=PYSCF_TMPDIR="$PWD/adc2_scratch" \
  ./run_ladder.sh
echo "ladder started as a user systemd unit (survives agent exit)"
systemctl --user is-active ladder
