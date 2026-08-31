#!/usr/bin/env bash
# Run the parse-time / behaviour guards.
#
#   .tests/guards/run.sh              every guard
#   .tests/guards/run.sh 36 37        only those numbers
#   .tests/guards/run.sh chimera      only guards whose filename matches
#   GUARD_KEEP_TMP=1 .tests/guards/run.sh 36    keep the workdir to inspect
#
# Each guard runs as its own process with its own temp directory, so one
# failing guard cannot corrupt another and any guard can be re-run alone.
#
# Portable to bash 3.2 (the /bin/bash macOS still ships): no mapfile, no
# associative arrays. Guard filenames never contain spaces, so plain
# word-splitting over a string list is safe here.
#
# `set -e` is deliberately NOT used. GitHub Actions injects a
# /usr/bin/clear_console call into the exit path of a `run:` block, and that
# helper exits non-zero on a non-TTY runner -- under -e it turned a fully
# passing run into a zero-output "exit code 1". Failures are tracked here
# explicitly instead.
set -uo pipefail

cd "$(dirname "$0")" || exit 1

ALL="$(ls -1 [0-9][0-9]_*.sh | sort)"

# Does guard file $1 match selector $2? A bare number selects that guard;
# anything else is a substring match on the filename.
guard_matches() {
  case "$2" in
    ''|*[!0-9]*)
      case "$1" in *"$2"*) return 0 ;; *) return 1 ;; esac ;;
    *)
      case "$1" in "$(printf '%02d_' "$2")"*) return 0 ;; *) return 1 ;; esac ;;
  esac
}

if [ "$#" -eq 0 ]; then
  SELECTED="$ALL"
else
  SELECTED=""
  for pat in "$@"; do
    for g in $ALL; do
      if guard_matches "$g" "$pat"; then
        case " $SELECTED " in
          *" $g "*) ;;
          *) SELECTED="$SELECTED $g" ;;
        esac
      fi
    done
  done
fi

SELECTED="$(echo $SELECTED)"
if [ -z "$SELECTED" ]; then
  echo "no guards matched: $*" >&2
  exit 2
fi

TOTAL=$(echo "$SELECTED" | wc -w | tr -d ' ')
FAILED=""
NFAILED=0
i=0
for g in $SELECTED; do
  i=$((i + 1))
  title="$(sed -n '2s/^# Guard [0-9]*: //p' "$g")"
  printf '[guards] %d/%d %s\n' "$i" "$TOTAL" "${title:-$g}"
  if ! "./$g"; then
    FAILED="$FAILED $g"
    NFAILED=$((NFAILED + 1))
  fi
done

echo
if [ "$NFAILED" -eq 0 ]; then
  echo "[guards] all $TOTAL passed"
  exit 0
fi
echo "[guards] $NFAILED of $TOTAL FAILED:"
for g in $FAILED; do
  echo "  $g"
  echo "      re-run: .tests/guards/$g"
done
exit 1
