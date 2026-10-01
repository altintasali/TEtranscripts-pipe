#!/usr/bin/env bash
# Guard 74: the splice-junctions screen is on by default
#
# Since 0.15.0 chimera.splice_junctions.enabled defaults to true: the screen
# only reads the SJ.out.tab files the main alignment already writes, so it
# costs no extra STAR pass. A config that never sets the key must schedule
# it; an explicit `enabled: false` must still turn it off.
#
# Run on its own:   .tests/guards/74_splice_junctions_screen_is_on_by_default.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# config/test.yaml sets the key explicitly; drop it to get the default
python3 - "$T" <<'PY'
import re, sys
T = sys.argv[1]
src = open("config/test.yaml").read()
block = re.search(r"^  splice_junctions:\n(?:(?:    .*|\s*)\n)*", src, re.M)
assert block, "no splice_junctions block in config/test.yaml"
body = block.group(0)
unset = re.sub(r"^    enabled: true\n", "", body, count=1, flags=re.M)
off = re.sub(r"^    enabled: true$", "    enabled: false", body, count=1, flags=re.M)
assert unset != body and off != body, "could not find splice_junctions.enabled"
open(f"{T}/unset.yaml", "w").write(src.replace(body, unset))
open(f"{T}/off.yaml", "w").write(src.replace(body, off))
PY

if ! snakemake --configfile "$T/unset.yaml" -n --cores 1 > "$T/unset.log" 2>&1; then
  echo "ERROR: dry run without splice_junctions.enabled failed"; tail -30 "$T/unset.log"; FAIL=1
elif ! grep -qE "^chimera_splice_junctions_classify\b" "$T/unset.log"; then
  echo "ERROR: with splice_junctions.enabled unset the screen must run by default"
  FAIL=1
fi

if ! snakemake --configfile "$T/off.yaml" -n --cores 1 > "$T/off.log" 2>&1; then
  echo "ERROR: dry run with splice_junctions.enabled: false failed"; tail -30 "$T/off.log"; FAIL=1
elif grep -qE "^chimera_splice_junctions_classify\b" "$T/off.log"; then
  echo "ERROR: an explicit enabled: false must still turn the screen off"
  FAIL=1
fi

exit $FAIL
