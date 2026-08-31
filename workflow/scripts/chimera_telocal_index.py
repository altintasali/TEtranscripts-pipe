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
import os
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


def telocal_te_id(key):
    """The individual TE insertion a TElocal key refers to, or None.

    This is the id te.bed and the chimera tables use as `te_id`, i.e. the TE
    GTF's transcript_id -- the FIRST colon field of a
    transcript_id:gene_id:family_id:class_id key. For a coordinate-style key
    (mghlab prebuilt annotations) the transcript_id is itself the
    `chrom:start:end(family:strand)` string, so that whole prefix is the id.
    """
    m = _LOCUS_RE.match(key)
    if m:
        close = key.find(")", m.end())
        return key[: close + 1] if close != -1 else None
    parts = key.split(":")
    return parts[0] if len(parts) >= 4 else None


def classify_telocal(key):
    """Return (family, class) from a TElocal TE key, or (None, None) for genes.

    Handles both key shapes TElocal produces (see TelocalIndex.build):
      chrom:start:end(family:strand):gene_id:family:class   (coordinate-style)
      transcript_id:gene_id:family_id:class_id              (the common case)
    In both, the last two colon-separated fields are family and class. Gene
    rows carry no colons and return (None, None).
    """
    m = _LOCUS_RE.match(key)
    if not m:
        # No coordinate prefix: the key is
        # transcript_id:gene_id:family_id:class_id (the common case -- see
        # TelocalIndex.build). The last two fields are still family and
        # class, so the same rule applies; anything with too few fields is a
        # gene row, which has no colons at all.
        parts = key.split(":")
        if len(parts) >= 4:
            return parts[-2], parts[-1]
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
    def build(cls, paths, open_read, locations=None):
        """Build from TElocal cntTables (paths), summing counts for the same
        locus key across all input files (samples).

        `locations` maps a cntTable row key to (chrom, start, end), read from
        telocal_locations.bed. It is the authoritative coordinate source
        because a TElocal key does NOT generally carry coordinates: the key
        is `transcript_id:gene_id:family_id:class_id`, and only when the TE
        GTF's transcript_id happens to be a coordinate string (the mghlab
        prebuilt annotations, e.g. `chr1:564318:564741(L1PA2:+)`) can they be
        parsed back out of it. An rmsk-derived GTF naming its insertions
        `L1PA2_dup1` -- which is the common case, and what our own
        build_telocal_index.py produces -- yields keys with no coordinates at
        all, so parse_telocal_locus returned None for every row, every row
        was skipped, and the index came out EMPTY. That made telocal_active
        "no" for every junction and the "TE locus is expressed" evidence tier
        silently unreachable. parse_telocal_locus is kept as a fallback for
        keys absent from the BED (a user-supplied .locInd built from a
        different annotation than ref.te_gtf).
        """
        locations = locations or {}
        raw = {}  # chrom -> {key: [start, end, key, family, cls, count]}
        n_rows = n_placed = 0
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
                    n_rows += 1
                    coords = locations.get(key) or parse_telocal_locus(key)
                    if coords is None:
                        continue
                    n_placed += 1
                    chrom, start, end = coords
                    family, cls_ = classify_telocal(key)
                    bucket = raw.setdefault(chrom, {})
                    row = bucket.get(key)
                    if row is None:
                        bucket[key] = [start, end, key, family, cls_, count]
                    else:
                        row[5] += count

        if n_rows and not n_placed:
            raise SystemExit(
                "error: none of the "
                f"{n_rows} TElocal rows could be given genomic coordinates, "
                "so the chimera TElocal index would be empty and every "
                "junction would be annotated telocal_active=no.\n"
                "This means the cntTable keys do not match any name in "
                "telocal_locations.bed and carry no coordinates themselves. "
                "The usual cause is a telocal.locind built from a DIFFERENT "
                "TE annotation than ref.te_gtf -- both must describe the "
                "same insertions, since the key "
                "(transcript_id:gene_id:family_id:class_id) is the join."
            )

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
        """Write atomically: full file to a temp name, then rename.

        A plain in-place write leaves a readable-but-TRUNCATED file if the
        job dies mid-write -- scheduler kill, quota, full filesystem. gzip
        and pickle have no length header, so nothing detects that until a
        downstream job calls pickle.load and gets a bare
        "EOFError: Ran out of input" from inside the stdlib, pointing at the
        reader rather than at the write that actually failed.
        os.replace is atomic within a filesystem, so the final path either
        does not exist or is a complete index -- never half of one.
        """
        tmp = f"{path}.tmp.{os.getpid()}"
        try:
            with gzip.open(tmp, "wb") as fh:
                pickle.dump(self, fh, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(tmp, path)
        except BaseException:
            if os.path.exists(tmp):
                os.remove(tmp)
            raise

    @classmethod
    def load(cls, path):
        try:
            with gzip.open(path, "rb") as fh:
                return pickle.load(fh)
        except (EOFError, pickle.UnpicklingError, gzip.BadGzipFile, OSError) as exc:
            size = os.path.getsize(path) if os.path.exists(path) else 0
            raise SystemExit(
                f"error: {path} is not a readable TElocal index "
                f"({type(exc).__name__}: {exc}); it is {size} bytes.\n"
                "The file is corrupt or truncated -- most likely a previous "
                "chimera_telocal_index job was killed or ran out of disk "
                "part-way through writing it.\n"
                "Delete it and re-run so it is rebuilt:\n"
                f"  rm {path}"
            ) from None

    def overlapping(self, chrom, start0, end0):
        """TE loci overlapping [start0, end0) on *chrom*.

        Correct for arbitrarily long/nested/overlapping loci, not just the
        common case -- bisects on the running max_end rather than start,
        matching classify_chimera_junctions.py's overlapping() rationale.
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
