#!/usr/bin/env bash
# Detach as a user systemd unit. setsid+nohup is NOT sufficient -- the agent
# harness kills its whole process tree on exit, and a setsid child is still a
# descendant. systemd reparents the job outside that tree, so it survives
# logout, the terminal closing, and the agent exiting.
#
# Attempt 4 (2 ranks x 24000 MB) was sized correctly -- both ranks stayed inside
# budget at 20.1 and 19.8 GB -- but was OOM-killed 4h19m in at IBatch 38/76
# because llama-server.service came up mid-run holding 20.5 GB, taking the box
# to 60.4 GB of 62. The kernel picked orca_mdci_mpi as the victim because a user
# service carries oom_score_adj=200. Hence the preflight guard below: this job
# needs the machine essentially to itself.
cd "$(dirname "$0")" || exit 1

NEED_GB=50
avail=$(free -g | awk 'NR==2{print $7}')
if [ "$avail" -lt "$NEED_GB" ]; then
  echo "REFUSING TO LAUNCH: ${avail} GB available, need >= ${NEED_GB} GB." >&2
  echo "Largest consumers:" >&2
  ps -eo rss,comm --sort=-rss | head -5 | awk 'NR>1{printf "  %7.1f GB  %s\n", $1/1048576, $2}' >&2
  echo "llama-server.service is the usual culprit: systemctl --user stop llama-server" >&2
  exit 1
fi
echo "preflight OK: ${avail} GB available"

systemctl --user reset-failed steom-tzvpd 2>/dev/null
# OOMScoreAdjust=0 asks systemd to undo the app.slice default of 200 so ORCA is
# not the preferred victim. Unprivileged users may not always be permitted to
# lower it; if systemd rejects the property the launch still proceeds without it.
systemd-run --user --unit=steom-tzvpd --working-directory="$PWD" \
  --property=OOMScoreAdjust=0 --property=ManagedOOMPreference=avoid \
  ./run_steom_tzvpd.sh \
  || systemd-run --user --unit=steom-tzvpd --working-directory="$PWD" ./run_steom_tzvpd.sh
echo "steom_def2_tzvpd running as user systemd unit 'steom-tzvpd'"
