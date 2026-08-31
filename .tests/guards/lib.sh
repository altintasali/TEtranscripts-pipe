#!/usr/bin/env bash
# Shared setup and fixtures for the guard scripts in this directory.
#
# Each guard is a standalone script: it gets its OWN temp directory, reports
# its own failures through $FAIL, and exits non-zero if any check failed. That
# is the whole point of the split -- previously all 49 lived in one 1,700-line
# inline block in ci.yml, so none could be run on its own and a failure meant
# scrolling a GitHub Actions log.
#
# Guards must not depend on each other's leftovers. Where several genuinely
# need the same fixture, it is a function here rather than an artifact the
# earlier guard happened to leave behind.

# Repo root, whatever directory the guard was invoked from.
GUARD_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Every guard runs from the repo root: the bodies reference config/test.yaml,
# .tests/resources/... and workflow/scripts/... as relative paths.
cd "$GUARD_REPO_ROOT" || exit 1

guard_init() {
  T="$(mktemp -d)"
  export T
  # shellcheck disable=SC2034  # read and incremented by the sourcing guard
  FAIL=0
  # Clean up unless the caller asked to keep the workdir for inspection.
  if [ -z "${GUARD_KEEP_TMP:-}" ]; then
    trap 'rm -rf "$T"' EXIT
  else
    trap 'echo "[guard] workdir kept: $T" >&2' EXIT
  fi
}

# --- shared fixtures --------------------------------------------------------

# Gzipped copies of the bundled test references, plus a config pointing at
# them. Used by every guard that exercises the .gz reference path.
fixture_gz_refs() {
  gzip -c .tests/resources/genome.fa > "$T/genome.fa.gz"
  gzip -c .tests/resources/genome.gtf > "$T/genome.gtf.gz"
  gzip -c .tests/resources/te_annotation.gtf > "$T/te_annotation.gtf.gz"
  sed -e "s|  fasta: .*|  fasta: $T/genome.fa.gz|" \
      -e "s|  gtf: .*|  gtf: $T/genome.gtf.gz|" \
      -e "s|  te_gtf: .*|  te_gtf: $T/te_annotation.gtf.gz|" \
      config/test.yaml > "$T/gz.yaml"
}

# The gene-TE evidence fixture and the merged table chimera_evidence.py builds
# from it. One pair per interesting case: called by both screens, canonical
# only, replicate-support only, assembly only, and a 999-read pair that must
# earn nothing (depth is the metric artifacts inflate most).
fixture_evidence() {
  mkdir -p "$T/ev"
  { printf 'event_id\tgene_id\tte_id\tte_subfamily\tte_family\tte_class\tcanonical\tchimera_type\ttelocal_active\tn_samples\ttotal_reads\n'
    printf 'j1\tGapdh\tL1PA2_dup1\tL1PA2\tL1\tLINE\tyes\tte_terminated\tyes\t6\t50\n'
    printf 'j2\tGapdh\tL1PA2_dup1\tL1PA2\tL1\tLINE\tno\tte_exonized\t.\t2\t10\n'
    printf 'j3\tActb\tAluY_dup9\tAluY\tAlu\tSINE\tno\tte_initiated\tyes\t3\t80\n'
    printf 'j4\tMyc\tL1MdA_dup4\tL1MdA\tL1\tLINE\tyes\tte_initiated\tno\t1\t5\n'
    printf 'j5\tTp53\tMIR_dup2\tMIR\tMIR\tSINE\tno\ttrans\tno\t1\t999\n'
  } | gzip -c > "$T/ev/j.tsv.gz"
  { printf 'transcript_id\tte_id\tte_subfamily\tte_family\tte_class\tmatched_gene_id\tstrand_match\tchimera_type\n'
    printf 'MSTRG.1.1\tL1PA2_dup1\tL1PA2\tL1\tLINE\tGapdh\tyes\tte_terminated\n'
    printf 'MSTRG.9.2\tSVA_dup3\tSVA\tSVA\tRetroposon\tSox2\tyes\tte_exonized\n'
  } | gzip -c > "$T/ev/a.tsv.gz"
  python3 workflow/scripts/chimera_evidence.py \
    --junction "$T/ev/j.tsv.gz" --assembly "$T/ev/a.tsv.gz" \
    --out "$T/ev/out.tsv.gz" > "$T/ev/fixture.log" 2>&1
}

# A directory that looks enough like a STAR index for the parse-time check in
# common.smk to accept it. Used by every guard that exercises
# star.build_index: false against a pre-existing index.
fixture_fake_star_index() {
  mkdir -p "$T/fake_star_index"
  touch "$T/fake_star_index/SA" "$T/fake_star_index/SAindex" \
        "$T/fake_star_index/Genome"
}
