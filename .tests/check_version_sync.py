#!/usr/bin/env python3
"""Assert workflow/environment.yaml and default-config/versions.yaml agree.

The same ~14 tools are pinned in two places:

  workflow/default-config/versions.yaml  what the workflow reports and what
                                         the generated per-tool conda envs use
  workflow/environment.yaml              the monolithic env people install

Both files carried a comment telling the reader to keep them in sync by hand.
That comment is the bug report: nothing enforced it, and a bump in one file
left the other silently stale -- so the report could claim STAR 2.7.11b while
the installed env ran something else.

This checks rather than generates on purpose. environment.yaml carries real
explanatory comments (why samtools is capped at 1.22 by STAR's htslib pin, why
TEtranscripts comes from PyPI and not bioconda) that are worth more where they
are than inside a template string. Generating the file would move them; this
keeps them and still makes drift impossible to merge.

Exit 0 when they agree, 1 with a diff when they do not.
"""
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
VERSIONS = ROOT / "workflow" / "default-config" / "versions.yaml"
ENVIRONMENT = ROOT / "workflow" / "environment.yaml"

# versions.yaml key -> package name as it appears in environment.yaml.
# Anything in versions.yaml with no entry here is reported as unmapped rather
# than silently skipped, so adding a tool cannot quietly bypass the check.
PACKAGE_NAMES = {
    "star": "star",
    "samtools": "samtools",
    "trim_galore": "trim-galore",
    "fastqc": "fastqc",
    "rseqc": "rseqc",
    "multiqc": "multiqc",
    "stringtie": "stringtie",
    "deseq2": "bioconductor-deseq2",
    "pheatmap": "r-pheatmap",
    "r_base": "r-base",
    "ucsc_gtftogenepred": "ucsc-gtftogenepred",
    "ucsc_genepredtobed": "ucsc-genepredtobed",
    "tetranscripts": "TEtranscripts",
    "telocal": "TElocal",
}


def env_pins(path):
    """package -> pinned version, across both conda and pip dependencies."""
    doc = yaml.safe_load(path.read_text())
    pins = {}
    for dep in doc.get("dependencies", []):
        if isinstance(dep, dict):
            for pip_dep in dep.get("pip", []) or []:
                m = re.match(r"^([A-Za-z0-9_.-]+)\s*==\s*(.+)$", pip_dep.strip())
                if m:
                    pins[m.group(1)] = m.group(2)
            continue
        m = re.match(r"^([A-Za-z0-9_.-]+)\s*=\s*([^=\s]+)$", str(dep).strip())
        if m:
            pins[m.group(1)] = m.group(2)
    return pins


def main():
    versions = yaml.safe_load(VERSIONS.read_text())["versions"]
    pins = env_pins(ENVIRONMENT)

    problems = []
    for key, want in sorted(versions.items()):
        package = PACKAGE_NAMES.get(key)
        if package is None:
            problems.append(
                f"  versions.yaml has '{key}' with no entry in PACKAGE_NAMES "
                f"({Path(__file__).name}) -- add one so it is checked"
            )
            continue
        got = pins.get(package)
        if got is None:
            problems.append(
                f"  {key}: versions.yaml pins {want!r}, but environment.yaml "
                f"has no pinned '{package}'"
            )
        elif str(got) != str(want):
            problems.append(
                f"  {key}: versions.yaml pins {want!r}, environment.yaml pins "
                f"{got!r} ({package})"
            )

    if problems:
        print("Tool versions disagree between the two files:\n")
        print("\n".join(problems))
        print(
            f"\nFix both:\n"
            f"  {VERSIONS.relative_to(ROOT)}\n"
            f"  {ENVIRONMENT.relative_to(ROOT)}"
        )
        return 1

    print(f"tool versions in sync ({len(versions)} pins checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
