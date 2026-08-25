"""Shared TElocal interval-index representation for the chimera<->TElocal
cross-reference (build_chimera_telocal_index.py / chimera_telocal_annotate.py).

Built ONCE from all samples' TElocal cntTables (build_chimera_telocal_index.py)
and reused by every per-sample chimera_telocal_annotate.py job, instead of
each of the N per-sample jobs re-parsing and re-summing all N samples'
cntTables from scratch (an Nx redundant rebuild -- the exact shape of the
bug that made telocal_counts OOM before its own two-pass rewrite, though
here the redundancy is across jobs rather than within one).

Per-chromosome storage uses parallel array.array('q'/'i', ...) columns for
the numeric fields (start/end/running-max-end/count) instead of a Python
dict-of-tuples: one boxed PyLong plus one dict entry per number is far
heavier than 8 raw bytes in a native array, and this index holds one row
per TE locus (potentially millions). family/class labels are deduplicated
into a small shared vocabulary (a handful of distinct TE families/classes
repeated across every locus) and stored as int indices rather than
repeating the same strings millions of times. Only the locus key string
itself (needed verbatim for the output's telocal_locus column) is
inherently unique per row and kept as a plain string list.
"""
import bisect
import gzip
import pickle
import re
from array import array as int_array

_LOCUS_RE = re.compile(r"^(.+):(\d+):(\d+)\(")


def parse_telocal_locus(key):
    """Extract (chrom, start, end) from a TElocal locus key.

    TElocal TE keys have the form
      chrom:start:end(family:strand):gene_id:family_id:class_id
    Gene keys (no parentheses in transcript portion) return None.
    """
    m = _LOCUS_RE.match(key)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def classify_telocal(key):
    """Return (family, class) from a TElocal TE key, or (None, None) for genes.

    TElocal TE keys have the form:
      chrom:start:end(family:strand):gene_id:family:class
    The last two colon-separated fields after the locus coordinates are
    family and class. Gene entries (no parentheses) return (None, None).
    """
    m = _LOCUS_RE.match(key)
    if not m:
        return None, None
    suffix = key[m.end():]
    # suffix is like "):L1PA2:L1:LINE/L1" or "):AluYb:SINE"
    # strip leading ')' or '):'
    suffix = suffix.lstrip("):").strip(":")
    parts = suffix.split(":")
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    if len(parts) == 1:
        return parts[0], None
    return None, None


class TelocalIndex:
    """Per-chromosome TE-locus interval index, columnar and vocab-deduped.

    overlapping() returns the same (start, end, key, family, cls, count)
    tuple shape the old dict-of-tuples version did, so callers (
    best_telocal_hit / annotate_event in chimera_telocal_annotate.py) are
    unchanged.
    """

    __slots__ = ("chroms", "vocab")

    def __init__(self):
        self.chroms = {}
        self.vocab = []

    @classmethod
    def build(cls, paths, open_read):
        """Build from TElocal cntTables (paths), summing counts for the same
        locus key across all input files (samples)."""
        raw = {}  # chrom -> {key: [start, end, key, family, cls, count]}
        for path in paths:
            with open_read(path) as fh:
                fh.readline()  # header: gene/TE \t <bam path>
                for line in fh:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    key = parts[0]
                    try:
                        count = int(parts[1])
                    except ValueError:
                        count = 0
                    coords = parse_telocal_locus(key)
                    if coords is None:
                        continue
                    chrom, start, end = coords
                    family, cls_ = classify_telocal(key)
                    bucket = raw.setdefault(chrom, {})
                    row = bucket.get(key)
                    if row is None:
                        bucket[key] = [start, end, key, family, cls_, count]
                    else:
                        row[5] += count

        index = cls()
        vocab_ids = {}

        def vocab_id(label):
            if label is None:
                return -1
            i = vocab_ids.get(label)
            if i is None:
                i = len(index.vocab)
                vocab_ids[label] = i
                index.vocab.append(label)
            return i

        for chrom, bucket in raw.items():
            feats = sorted(bucket.values(), key=lambda r: (r[0], r[1]))
            n = len(feats)
            starts = int_array("q", (r[0] for r in feats))
            ends = int_array("q", (r[1] for r in feats))
            counts = int_array("q", (r[5] for r in feats))
            fam_ids = int_array("i", (vocab_id(r[3]) for r in feats))
            cls_ids = int_array("i", (vocab_id(r[4]) for r in feats))
            keys = [r[2] for r in feats]
            max_end = int_array("q", bytes(8 * n))
            running = -1
            for i, e in enumerate(ends):
                if e > running:
                    running = e
                max_end[i] = running
            index.chroms[chrom] = {
                "start": starts,
                "end": ends,
                "max_end": max_end,
                "count": counts,
                "family": fam_ids,
                "class": cls_ids,
                "key": keys,
            }
        return index

    def save(self, path):
        with gzip.open(path, "wb") as fh:
            pickle.dump(self, fh, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, path):
        with gzip.open(path, "rb") as fh:
            return pickle.load(fh)

    def overlapping(self, chrom, start0, end0):
        """TE loci overlapping [start0, end0) on *chrom*.

        Correct for arbitrarily long/nested/overlapping loci, not just the
        common case -- bisects on the running max_end rather than start,
        matching parse_chimeric_junctions.py's overlapping() rationale.
        """
        cols = self.chroms.get(chrom)
        if cols is None:
            return []
        starts, ends, max_end, counts = cols["start"], cols["end"], cols["max_end"], cols["count"]
        fam_ids, cls_ids, keys = cols["family"], cols["class"], cols["key"]
        lo = bisect.bisect_right(max_end, start0)
        hits = []
        for i in range(lo, len(keys)):
            if starts[i] >= end0:
                break
            if ends[i] > start0:
                fam = self.vocab[fam_ids[i]] if fam_ids[i] >= 0 else None
                cls_ = self.vocab[cls_ids[i]] if cls_ids[i] >= 0 else None
                hits.append((starts[i], ends[i], keys[i], fam, cls_, counts[i]))
        return hits
