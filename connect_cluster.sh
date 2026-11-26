#!/usr/bin/env bash
# connect_cluster.sh - open Fluorite Ordering Studio on your cluster files from your laptop.
#
# Run this on your LAPTOP. It logs in to the cluster once (password / 2FA asked once), starts the app on the
# login node it landed on (or reuses the one already running there), forwards a local port to it through the
# same SSH connection, and prints (and opens) a http://127.0.0.1 link. The app, the files and SLURM stay on
# the cluster, so every page works, including Submit. Ctrl+C disconnects and stops the app on the cluster
# (add -k to leave it running). Jobs started from the app keep running either way.
#
#   ./connect_cluster.sh you@login.cluster.org
#   ./connect_cluster.sh -a ~/fluorite-ordering-studio -d /scratch/me/runs -s 'module load Python' you@host
#   ./connect_cluster.sh --stop you@host          # stop the app on that login node
#
# Options (all can be saved in ~/.fluorite_studio/cluster.conf as HOST= APP= DIR= SETUP= PY= PORT= KEEP=):
#   -a DIR    app folder on the cluster, where start_studio.py is   (default ~/fluorite-ordering-studio)
#   -d DIR    folder the app opens in, e.g. your run folder         (default: your home on the cluster)
#   -s CMD    shell commands run before python on the cluster, e.g. 'module load foss Python'
#   -P PY     python on the cluster                                 (default python3)
#   -p PORT   first local port to try on this laptop                (default 8765)
#   -k        keep the app running on the cluster after disconnecting
#   --stop    stop the app on the cluster and exit
#   -o OPT    extra ssh option, e.g. -o 'Port=2222' or -o 'IdentityFile=~/.ssh/id_cluster' (repeatable)
#   -n        do not open a browser
set -u

HOST="" APP="~/fluorite-ordering-studio" DIR="" SETUP="" PY="python3" PORT=8765 KEEP=0 NOOPEN=0 STOP=0
CONF="${FO_STUDIO_CLUSTER_CONF:-$HOME/.fluorite_studio/cluster.conf}"
# shellcheck disable=SC1090
[ -f "$CONF" ] && . "$CONF"
SSHOPTS=()
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit "${1:-0}"; }
while [ $# -gt 0 ]; do
  case "$1" in
    -a) APP="$2"; shift 2;;        -d) DIR="$2"; shift 2;;
    -s) SETUP="$2"; shift 2;;      -P) PY="$2"; shift 2;;
    -p) PORT="$2"; shift 2;;       -o) SSHOPTS+=(-o "$2"); shift 2;;
    -k) KEEP=1; shift;;            -n) NOOPEN=1; shift;;
    --stop) STOP=1; shift;;        -h|--help) usage 0;;
    -*) echo "unknown option $1" >&2; usage 2;;
    *) HOST="$1"; shift;;
  esac
done
[ -n "$HOST" ] || { echo "give the cluster as user@host (or HOST= in $CONF)" >&2; usage 2; }
case "$PORT" in ''|*[!0-9]*) echo "-p needs a number" >&2; exit 2;; esac

CTL_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fo_studio.XXXXXX") || exit 1
CTL="$CTL_DIR/ctl"
SSH=(ssh "${SSHOPTS[@]}" -o ControlPath="$CTL")
q() { printf '%q' "$1"; }

# 1. one login; every later command and the port forward go through it, so they reach the same login node
echo "connecting to $HOST ..."
if ! ssh "${SSHOPTS[@]}" -o ControlMaster=yes -o ControlPath="$CTL" -o ServerAliveInterval=30 -fN "$HOST"; then
  rmdir "$CTL_DIR" 2>/dev/null; echo "ssh login to $HOST failed" >&2; exit 1
fi
close_master() { "${SSH[@]}" -O exit "$HOST" >/dev/null 2>&1; rm -rf "$CTL_DIR"; }

# 2. start the app on that login node (or find the one already running there)
REMOTE=$(cat <<'EOS'
set -u
ACT=$1 APP=$2 DIR=$3 SETUP=$4 PY=$5
case "$APP" in "~"|"~/"*) APP="$HOME${APP#\~}";; esac
case "$DIR" in "~"|"~/"*) DIR="$HOME${DIR#\~}";; esac
ST="$HOME/.fluorite_studio"; mkdir -p "$ST"; chmod 700 "$ST"
NODE=$(hostname); PIDF="$ST/server_$NODE.pid"; LOG="$ST/server_$NODE.log"
alive() { [ -f "$PIDF" ] && kill -0 "$(cat "$PIDF")" 2>/dev/null; }
if [ "$ACT" = stop ]; then
  if alive; then kill "$(cat "$PIDF")"; echo "stopped the app on $NODE"; else echo "no app running on $NODE"; fi
  rm -f "$PIDF"; exit 0
fi
if alive && grep -q '^  open    : ' "$LOG" 2>/dev/null; then
  echo "REUSED"
else
  [ -f "$APP/start_studio.py" ] || { echo "ERROR no start_studio.py in $APP on $NODE (use -a)"; exit 1; }
  [ -z "$DIR" ] || [ -d "$DIR" ] || { echo "ERROR folder $DIR not found on $NODE (use -d)"; exit 1; }
  if [ -n "$SETUP" ]; then eval "$SETUP" || { echo "ERROR setup command failed: $SETUP"; exit 1; }; fi
  command -v "$PY" >/dev/null 2>&1 || { echo "ERROR $PY not found on $NODE (use -P or -s 'module load ...')"; exit 1; }
  cd "$APP" || exit 1
  set -- "$PY" -u start_studio.py; [ -n "$DIR" ] && set -- "$@" --dir "$DIR"
  if command -v setsid >/dev/null 2>&1; then setsid nohup "$@" >"$LOG" 2>&1 </dev/null &
  else nohup "$@" >"$LOG" 2>&1 </dev/null & fi
  echo $! >"$PIDF"
  for _ in $(seq 1 60); do
    grep -q '^  open    : ' "$LOG" 2>/dev/null && break
    kill -0 "$(cat "$PIDF")" 2>/dev/null || break
    sleep 0.5
  done
  grep -q '^  open    : ' "$LOG" 2>/dev/null || { echo "ERROR the app did not start on $NODE:"; tail -n 20 "$LOG"; rm -f "$PIDF"; exit 1; }
  echo "STARTED"
fi
echo "NODE $NODE"
echo "PID $(cat "$PIDF")"
sed -n 's/^  open    : /URL /p' "$LOG" | tail -n 1
EOS
)
run_remote() {   # $1 = start|stop
  "${SSH[@]}" "$HOST" "bash -l -c $(q "$REMOTE") fo_studio $(q "$1") $(q "$APP") $(q "$DIR") $(q "$SETUP") $(q "$PY")"
}

if [ "$STOP" = 1 ]; then run_remote stop; rc=$?; close_master; exit $rc; fi

OUT=$(run_remote start); rc=$?
URL=$(printf '%s\n' "$OUT" | sed -n 's/^URL //p' | tail -n 1)
if [ $rc -ne 0 ] || [ -z "$URL" ]; then
  printf '%s\n' "$OUT" | grep -v '^\(NODE\|PID\|URL\) ' >&2; close_master; exit 1
fi
NODE=$(printf '%s\n' "$OUT" | sed -n 's/^NODE //p')
STATE=$(printf '%s\n' "$OUT" | grep -m1 -E '^(STARTED|REUSED)$')
RPORT=$(printf '%s\n' "$URL" | sed -n 's|^http://[^/]*:\([0-9]*\)/.*|\1|p')
TOKEN=$(printf '%s\n' "$URL" | sed -n 's|.*[?&]token=\([^&]*\).*|\1|p')
[ -n "$RPORT" ] && [ -n "$TOKEN" ] || { echo "could not read the app link: $URL" >&2; close_master; exit 1; }

# 3. forward the first free local port to the app's port on that node
LPORT=""
for p in $(seq "$PORT" $((PORT + 19))); do
  if "${SSH[@]}" -o ExitOnForwardFailure=yes -O forward -L "127.0.0.1:$p:127.0.0.1:$RPORT" "$HOST" >/dev/null 2>&1; then
    LPORT=$p; break
  fi
done
[ -n "$LPORT" ] || { echo "no free local port in $PORT-$((PORT + 19))" >&2; close_master; exit 1; }
LINK="http://127.0.0.1:$LPORT/?token=$TOKEN"

cleanup() {
  trap - INT TERM HUP
  echo
  if [ "$KEEP" = 1 ]; then echo "disconnected; the app keeps running on $NODE (stop it with --stop)"
  else run_remote stop; fi
  close_master; exit 0
}
trap cleanup INT TERM HUP

echo "Fluorite Ordering Studio on $NODE ($(printf %s "$STATE" | tr A-Z a-z), cluster port $RPORT)"
echo "  open    : $LINK"
echo "  stop    : Ctrl+C$( [ "$KEEP" = 1 ] && echo ' (the app keeps running, -k)' || echo ' (stops the app; running jobs keep running)')"
if [ "$NOOPEN" = 0 ]; then
  if command -v xdg-open >/dev/null 2>&1; then xdg-open "$LINK" >/dev/null 2>&1 &
  elif command -v open >/dev/null 2>&1; then open "$LINK" >/dev/null 2>&1 & fi
fi

# 4. stay connected until Ctrl+C or the connection drops
while "${SSH[@]}" -O check "$HOST" >/dev/null 2>&1; do sleep 2 & wait $!; done
echo "connection to $HOST lost; run the same command again to reconnect (the app is still on $NODE)"
rm -rf "$CTL_DIR"; exit 1
