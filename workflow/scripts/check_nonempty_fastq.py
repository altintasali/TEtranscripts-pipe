"""Fail if a trimmed FASTQ file is empty (zero reads).

Called from trim_galore's shell to catch the case where TrimGalore
discards all reads, producing a valid-but-empty gzip file.  The
pipeline would otherwise silently propagate empty data through STAR,
TEcount, and every downstream stage.

Usage: check_nonempty_fastq.py FILE1 [FILE2 ...]
Non-existent paths are silently skipped (only PE or SE outputs exist).
"""
import gzip
import os
import sys


def has_reads(path):
    """Return True if the FASTQ/GZ file contains at least one read."""
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("@"):
                    return True
    except Exception:
        pass
    return False


def main():
    existing = [p for p in sys.argv[1:] if os.path.exists(p)]
    if not existing:
        print(
            "ERROR: trimming produced no output files", file=sys.stderr
        )
        sys.exit(1)
    for path in existing:
        if has_reads(path):
            print(f"OK: {path} has reads", file=sys.stderr)
            return
    # None of the existing files had reads
    print(
        f"ERROR: trimming produced empty FASTQ(s) "
        f"(all reads removed): {', '.join(existing)}",
        file=sys.stderr,
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
