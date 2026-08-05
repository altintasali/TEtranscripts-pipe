#!/usr/bin/env python3
"""Aggregate Snakemake benchmark files into a self-describing MultiQC custom
content file rendered as the "Resource usage" section of the workflow's
MultiQC report.

Each input is a Snakemake benchmark .txt (tab-separated rows with the header
`s  seconds  threads  cpu_percent  max_rss  bytes_read  bytes_written`,
max_rss in KB). Files are grouped by rule (the benchmark subdirectory name),
and per-rule summary stats -- job count, mean/max wall time, mean/max peak
RSS, and mean CPU load -- are written as a table under the fixed id
"resource_usage" so the module_order in multiqc_config.yaml places it after
the RSeQC section.
"""
import json
import os
import statistics
from collections import defaultdict

try:
    # Only defined when run by Snakemake's `script:` directive.
    snakemake
except NameError:  # pragma: no cover
    snakemake = None


def _parse_benchmark(path):
    """Yield one dict per run row in a Snakemake benchmark file."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if header is None:
                header = line.split("\t")
                continue
            yield dict(zip(header, line.split("\t")))


def main():
    rows_by_rule = defaultdict(list)
    for path in snakemake.input:
        rule = os.path.basename(os.path.dirname(path))
        rows_by_rule[rule].extend(_parse_benchmark(path))

    data = {}
    for rule in sorted(rows_by_rule):
        rows = rows_by_rule[rule]
        walltimes = [float(r.get("s") or 0) for r in rows]
        rss_kb = [float(r.get("max_rss") or 0) for r in rows]
        cpu = [float(r.get("cpu_percent") or 0) for r in rows]
        data[rule] = {
            "n": len(rows),
            "walltime_mean_h": round(statistics.mean(walltimes) / 3600.0, 3),
            "walltime_max_h": round(max(walltimes) / 3600.0, 3),
            "rss_mean_gb": round(statistics.mean(rss_kb) / (1024.0 * 1024.0), 3),
            "rss_max_gb": round(max(rss_kb) / (1024.0 * 1024.0), 3),
            "mean_load": round(statistics.mean(cpu) / 100.0, 3),
        }

    summary = {
        "id": "resource_usage",
        "section_name": "Resource usage",
        "description": (
            "Per-rule job count, wall time, peak memory and CPU load, "
            "aggregated from Snakemake's benchmark files "
            "(results/pipeline_info/benchmarks/). Useful for sizing "
            "resources on your cluster before a full run."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "resource_usage_table",
            "title": "Resource usage",
            "col1_header": "Rule",
            "sort_rows": False,
        },
        "headers": {
            "n": {
                "title": "N",
                "description": "Number of jobs (samples or files)",
                "format": "{:,d}",
                "min": 0,
            },
            "walltime_mean_h": {
                "title": "Wall time mean (h)",
                "format": "{:.3f}",
                "min": 0,
            },
            "walltime_max_h": {
                "title": "Wall time max (h)",
                "format": "{:.3f}",
                "min": 0,
            },
            "rss_mean_gb": {
                "title": "Peak RSS mean (GB)",
                "format": "{:.3f}",
                "min": 0,
            },
            "rss_max_gb": {
                "title": "Peak RSS max (GB)",
                "format": "{:.3f}",
                "min": 0,
            },
            "mean_load": {
                "title": "Mean CPU load",
                "description": "Mean cpu_percent / 100 across jobs",
                "format": "{:.3f}",
                "min": 0,
            },
        },
        "data": data,
    }

    with open(snakemake.output[0], "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
