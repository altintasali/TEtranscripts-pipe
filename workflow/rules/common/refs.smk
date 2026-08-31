# Reference files (incl. gzipped), sjdboverhang, star index location.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

# -----------------------------------------------------------------------------
# Reference files: transparently support gzipped fasta/gtf/te_gtf.
# STAR's --genomeFastaFiles/--sjdbGTFfile and TEcount/TEtranscripts'
# --GTF/--TE all expect plain-text files, so any ref.* path ending in .gz is
# decompressed once via the gunzip_reference rule (ref.smk); everything else
# in the workflow refers to the resolved (always-uncompressed) path below,
# never to config["ref"][...] directly.
#
# Decompressed files land in ref.decompressed_dir. The default is a
# directory under the results tree (results/pipeline_info/ref_decompressed):
# unlike /tmp, that lives on the same shared filesystem the snakemake
# scheduler and every cluster job can see, so gunzip_reference's output is
# always visible -- a node-local /tmp path silently breaks cluster runs with
# "output ... missing locally, parent dir not present" (the file is written
# on the compute node, but the scheduler checks for it on the submission
# node). The files are temp()-wrapped (see the gunzip_reference rule) so they
# are deleted after downstream rules consume them. Set ref.decompressed_dir
# to a scratch dir on shared storage if you would rather keep them.
#
# Because a node-local dir (or an empty value, which would resolve to
# "/<stem>" at the filesystem root) can only break runs that actually
# decompress something, we hard-fail at parse time when gzipped references
# are present -- before any job is submitted -- instead of letting the first
# gunzip_reference job die hours into a queue. Configs with plain refs never
# reach gunzip_reference, so an explicit decompressed_dir is dead config
# there and left alone.
# -----------------------------------------------------------------------------
DECOMPRESS_DIR = config["ref"].get(
    "decompressed_dir", "results/pipeline_info/ref_decompressed"
)

_HAS_GZ_REFS = any(
    str(config["ref"].get(k, "")).endswith(".gz")
    for k in ("fasta", "gtf", "te_gtf")
)
if "decompressed_dir" in config.get("ref", {}) and _HAS_GZ_REFS:
    _local_tmp = tempfile.gettempdir().rstrip(os.sep)
    _node_local = (
        DECOMPRESS_DIR == "/tmp"
        or DECOMPRESS_DIR.startswith("/tmp/")
        or DECOMPRESS_DIR == "/var/tmp"
        or DECOMPRESS_DIR.startswith("/var/tmp/")
        or "$TMPDIR" in DECOMPRESS_DIR
        or (_local_tmp and DECOMPRESS_DIR.startswith(_local_tmp))
    )
    if _node_local or not DECOMPRESS_DIR.strip():
        raise WorkflowError(
            f"Gzipped reference files are used, but ref.decompressed_dir is "
            f"{DECOMPRESS_DIR!r}{' (empty)' if not DECOMPRESS_DIR.strip() else ''} "
            f"-- a node-local temp directory. gunzip_reference writes its "
            f"output there, and on a multi-node cluster each job runs on a "
            f"compute node with its own /tmp that the snakemake scheduler "
            f"cannot see ('MissingOutputException ... (missing locally, parent "
            f"dir not present)'). Delete the ref.decompressed_dir line to use "
            f"the default results/pipeline_info/ref_decompressed (shared "
            f"storage, temp()-cleaned), or point it at a scratch directory on "
            f"shared storage."
        )

REFERENCE_GZ_SOURCES = {}  # decompressed stem -> original .gz path, for gunzip_reference


def _resolve_ref_path(key):
    path = config["ref"][key]
    if str(path).endswith(".gz"):
        stem = os.path.basename(str(path))[:-3]  # strip trailing ".gz"
        decompressed = f"{DECOMPRESS_DIR}/{stem}"
        if stem in REFERENCE_GZ_SOURCES and REFERENCE_GZ_SOURCES[stem] != path:
            raise ValueError(
                f"Two different ref.* .gz files decompress to the same "
                f"filename '{stem}' ({decompressed}) -- "
                "rename one of them so their basenames are unique."
            )
        REFERENCE_GZ_SOURCES[stem] = path
        return decompressed
    return path


# FASTA is only needed when building the index; when build_index is false
# the user provides a pre-built index and ref.fasta is irrelevant.
_STAR_BUILD_EARLY = config["star"].get("build_index", True)
if _STAR_BUILD_EARLY:
    FASTA = _resolve_ref_path("fasta")
else:
    _raw_fasta = (config["ref"].get("fasta") or "").strip()
    if _raw_fasta:
        import warnings as _warn
        _warn.warn(
            f"ref.fasta ({_raw_fasta}) is set but star.build_index is false "
            f"-- the FASTA file will not be used. Remove it from config.yaml "
            f"to suppress this warning.",
            stacklevel=2,
        )
    FASTA = None  # not needed; the external index is used as-is

GTF = _resolve_ref_path("gtf")
TE_GTF = _resolve_ref_path("te_gtf")

# References are required for everything downstream (index, gene model,
# TEcount), so fail fast at parse time on a missing/typo'd path instead of
# letting the first rule that touches it die hours into a queue. For .gz refs
# we check the .gz source: the decompressed target is a temp() rule output
# that legitimately doesn't exist until gunzip_reference runs.
_missing_refs = []
for _key, _resolved in (("gtf", GTF), ("te_gtf", TE_GTF)):
    _source = config["ref"][_key]
    _check = _source if str(_source).endswith(".gz") else _resolved
    if not os.path.exists(_check):
        _missing_refs.append(
            f"ref.{_key}: {_check} (exists: False, resolves to: "
            f"{os.path.abspath(_check)})"
        )
if FASTA is not None:
    _fasta_src = config["ref"]["fasta"]
    _fasta_check = _fasta_src if str(_fasta_src).endswith(".gz") else FASTA
    if not os.path.exists(_fasta_check):
        _missing_refs.append(
            f"ref.fasta: {_fasta_check} (exists: False, resolves to: "
            f"{os.path.abspath(_fasta_check)})"
        )
if _missing_refs:
    raise WorkflowError(
        "Missing reference file(s) from config.yaml:\n  "
        + "\n  ".join(_missing_refs)
        + "\nFix the paths in config.yaml (they're resolved relative to the "
        "directory you run snakemake from, usually the repo root)."
    )

# Guard against empty GTF files: a file that exists but contains only
# comments or blank lines would silently produce empty BED tracks and
# zero-row count matrices, leading to broken downstream analysis.
def _count_gtf_features(path):
    """Return the number of non-comment, non-empty lines in a GTF file."""
    import gzip as _gzip

    _open = _gzip.open if str(path).endswith(".gz") else open
    _count = 0
    try:
        with _open(path, "rt") as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    _count += 1
    except Exception:
        return -1  # unreadable — let the file-exists check catch it
    return _count


_empty_gtf_warnings = []
for _key, _resolved in (("gtf", GTF), ("te_gtf", TE_GTF)):
    _source = config["ref"][_key]
    _check = _source if str(_source).endswith(".gz") else _resolved
    if os.path.exists(_check):
        _n = _count_gtf_features(_check)
        if _n == 0:
            _empty_gtf_warnings.append(
                f"ref.{_key}: {_check} exists but contains no feature lines "
                f"(empty GTF). Downstream results will be empty."
            )
        elif _n < 0:
            _empty_gtf_warnings.append(
                f"ref.{_key}: {_check} exists but could not be read to "
                f"validate content."
            )
if _empty_gtf_warnings:
    logger.warning(
        "Empty or unreadable GTF file(s):\n  "
        + "\n  ".join(_empty_gtf_warnings)
        + "\nDownstream BED tracks, count matrices, and QC views will be "
        "empty for these annotations."
    )

# -----------------------------------------------------------------------------
# STAR sjdbOverhang: auto-detected from the sample sheet's fastq files
# (max read length - 1, STAR's own recommendation) unless the user pins an
# explicit integer in config.yaml. Detection peeks at the first few reads of
# every fastq_1/fastq_2 file, gzip-aware, so it works whether reads are
# compressed or not.
# -----------------------------------------------------------------------------
def _open_maybe_gz(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _fastq_read_length(path, n_reads=20):
    """Max read length in the first n_reads of a fastq(.gz) file, or 0 if
    the file doesn't exist / can't be read (e.g. not downloaded yet)."""
    if not path or not os.path.exists(path):
        return 0
    max_len = 0
    try:
        with _open_maybe_gz(path) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:  # sequence line of each fastq record
                    max_len = max(max_len, len(line.strip()))
                if i >= n_reads * 4:
                    break
    except OSError:
        return 0
    return max_len


def _resolve_sjdb_overhang():
    configured = config["ref"].get("sjdb_overhang", "auto")
    if str(configured).strip().lower() != "auto":
        return int(configured), "config.yaml (explicit)"
    paths = []
    for s in SAMPLES:
        for read in (1, 2):
            paths.extend(sample_fastqs(s, read))
    max_len = max((_fastq_read_length(p) for p in paths), default=0)
    if max_len == 0:
        cwd = os.getcwd()
        checked = "\n".join(
            f"  - {p}  (exists: {os.path.exists(p)}, resolves to: "
            f"{os.path.abspath(p)})"
            for p in paths
        ) or "  (no fastq_1/fastq_2 values found in the sample sheet at all)"
        raise ValueError(
            "ref.sjdb_overhang is 'auto' but no readable fastq file was "
            f"found to auto-detect read length from. Current working "
            f"directory (paths are resolved relative to this): {cwd}\n"
            f"Paths checked from {config['samples']}:\n{checked}\n"
            "Either run Snakemake from the directory these paths are "
            "relative to (usually the repo root), fix the paths in the "
            "sample sheet, or set ref.sjdb_overhang to an explicit integer "
            "in config.yaml to skip auto-detection entirely."
        )
    return max_len - 1, f"auto-detected (max read length {max_len} - 1)"


SJDB_OVERHANG, _SJDB_OVERHANG_SOURCE = _resolve_sjdb_overhang()

# Record how key values were resolved -- a lightweight, always-current log of
# workflow parse-time decisions, separate from the per-rule logs in
# results/pipeline_info/logs/.
os.makedirs("results/pipeline_info/logs", exist_ok=True)
with open("results/pipeline_info/logs/config_resolution.log", "w") as fh:
    fh.write(f"sjdb_overhang = {SJDB_OVERHANG} ({_SJDB_OVERHANG_SOURCE})\n")
    fh.write(f"decompressed reference directory = {DECOMPRESS_DIR}\n")
    if REFERENCE_GZ_SOURCES:
        fh.write("gzipped reference files to be decompressed:\n")
        for stem, gz_path in REFERENCE_GZ_SOURCES.items():
            fh.write(f"  {gz_path} -> {DECOMPRESS_DIR}/{stem}\n")
    else:
        fh.write("no gzipped reference files (fasta/gtf/te_gtf) detected\n")

# -----------------------------------------------------------------------------
# Where the STAR index lives.
#
# build_index=true  (default): the index is built (or reused if already
#   present) at the path given by star.index, or results/star_index if unset.
#
# build_index=false: the user provides a pre-built external index via
#   star.index and it is used in place, as a plain rule INPUT.  No rule
#   declares it as an output, so Snakemake never owns -- and never deletes
#   -- it.  See the longer note at the else: branch below.
# -----------------------------------------------------------------------------
_STAR_INDEX_RAW = (config["star"].get("index") or "").strip()
STAR_BUILD_INDEX = config["star"].get("build_index", True)
