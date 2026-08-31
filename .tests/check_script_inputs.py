#!/usr/bin/env python3
"""Assert every rule that runs a script declares that script as an input.

Snakemake's `code` rerun-trigger hashes a rule's shell command STRING, not the
file the command names. So a rule written as

    shell: "python3 {SCRIPTS_DIR}/foo.py ..."

is NOT re-run when foo.py's contents change: snakemake reports "Nothing to be
done" and the stale output survives. Measured -- editing a script's behaviour
left the previous output in place with no warning.

Declaring the script as an input fixes it (the reason becomes "Updated input
files: foo.py"), so this checks that none of them regress. The `script:`
directive is already tracked by snakemake and needs no declaration.

Exit 0 when every rule declares its script, 1 with the offenders listed.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RULE = re.compile(r"^(?P<indent> *)rule (?P<name>[a-z_0-9]+):\n", re.M)


def main():
    problems = []
    checked = 0
    for f in sorted((ROOT / "workflow" / "rules").rglob("*.smk")):
        src = f.read_text()
        starts = [(m.start(), m.group("name")) for m in RULE.finditer(src)]
        for i, (start, name) in enumerate(starts):
            end = starts[i + 1][0] if i + 1 < len(starts) else len(src)
            body = src[start:end]
            invoked = set(re.findall(r"\{SCRIPTS_DIR\}/([A-Za-z0-9_.-]+)", body))
            if not invoked:
                continue
            checked += 1
            declared = set(re.findall(
                r'script\s*=\s*f?"\{SCRIPTS_DIR\}/([A-Za-z0-9_.-]+)"', body))
            # a rule may legitimately name the script only in its declaration
            undeclared = invoked - declared
            if undeclared:
                problems.append(
                    f"  {f.name}:{name} runs {sorted(undeclared)} without "
                    f"declaring it as an input")

    if problems:
        print("Rules that run a script without declaring it:\n")
        print("\n".join(problems))
        print("\nAdd it to the rule's input block and use {input.script}:\n"
              '    input:\n'
              '        script=f"{SCRIPTS_DIR}/<name>",\n'
              "Otherwise editing that script leaves stale outputs in place.")
        return 1

    print(f"every script-running rule declares its script ({checked} checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
