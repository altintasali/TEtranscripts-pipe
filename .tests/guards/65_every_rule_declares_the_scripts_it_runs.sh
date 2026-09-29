#!/usr/bin/env bash
# Guard 65: every rule declares the workflow/scripts/ file it runs as an input
#
# Snakemake does not reliably notice when only a script changes:
#   - for shell: rules, the "code" rerun trigger hashes the shell command
#     STRING, not the file it names;
#   - for script: rules, it is not something this pipeline can count on
#     across Snakemake versions either.
# Declaring the script as an input makes the ordinary mtime trigger rerun the
# rule (and everything downstream) whenever the script is edited. Without it,
# a fix to a report script leaves the old output -- and the old report -- in
# place silently. That happened on a real run: three committed fixes never
# reached the rebuilt MultiQC report.
#
# The check is static: every rule/checkpoint block (nested ones included) is
# scanned, comments stripped, and every workflow/scripts/ file it names --
# via a script: directive or literally anywhere in the rule -- must also be
# named in that rule's input: section. Files that are not in workflow/scripts/
# (e.g. RSeQC's own read_distribution.py) are ignored.
#
# Run on its own:   .tests/guards/65_every_rule_declares_the_scripts_it_runs.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

python3 - <<'PY' || FAIL=1
import glob
import os
import re
import sys

scripts = {os.path.basename(p) for p in glob.glob("workflow/scripts/*")
           if os.path.isfile(p)}
head = re.compile(r"^(\s*)(?:rule|checkpoint)\s+(\w+)\s*:")
directive = re.compile(r"^\s*(\w+):\s*(.*)$")
name_re = re.compile(r"[\w.-]+\.(?:py|R|sh)\b")


def strip_comment(line):
    # good enough for .smk: no rule here puts a '#' inside a string literal
    # that names a script
    return line.split("#", 1)[0]


problems = []
n_rules = 0
for path in sorted(glob.glob("workflow/rules/**/*.smk", recursive=True)):
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        m = head.match(lines[i])
        if not m:
            i += 1
            continue
        indent, rule = len(m.group(1)), m.group(2)
        body = []
        j = i + 1
        while j < len(lines):
            ln = lines[j]
            if ln.strip() and (len(ln) - len(ln.lstrip())) <= indent:
                break
            body.append(strip_comment(ln))
            j += 1
        n_rules += 1

        # split the body into directive sections (input:, output:, shell:, ...)
        sections, current = {}, None
        dir_indent = None
        for ln in body:
            if not ln.strip():
                continue
            ind = len(ln) - len(ln.lstrip())
            d = directive.match(ln)
            if d and (dir_indent is None or ind <= dir_indent):
                dir_indent = ind
                current = d.group(1)
                sections.setdefault(current, []).append(d.group(2))
            elif current:
                sections[current].append(ln)
        text = {k: "\n".join(v) for k, v in sections.items()}
        declared = set(name_re.findall(text.get("input", "")))
        used = {n for n in name_re.findall("\n".join(body))
                if n in scripts}
        missing = sorted(used - declared)
        if missing:
            problems.append(f"{path}: rule {rule} runs {', '.join(missing)} "
                            "but does not declare it as an input")
        i = j

if n_rules < 50:
    print(f"ERROR: only {n_rules} rules found -- the scanner is broken")
    sys.exit(1)
for p in problems:
    print("ERROR:", p)
sys.exit(1 if problems else 0)
PY

exit $FAIL
