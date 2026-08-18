"""Gzip-aware file helpers for TEtranscripts-pipe scripts.

Pipeline text outputs use .gz extensions to save disk space.  These helpers
let every script transparently read and write compressed files without
scattering ``gzip.open`` conditionals throughout the codebase.
"""
import gzip


def open_read(path):
    """Open *path* for reading, transparently decompressing if it ends in .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def open_write(path):
    """Open *path* for writing, transparently compressing if it ends in .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")
