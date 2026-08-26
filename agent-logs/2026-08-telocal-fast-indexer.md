# TElocal fast indexer — session log

Branch: `fix/fast-telocal-indexer` (based on `feature/telocal`, deliberately
*not* on `feature/chimera-assembly` — keeps this fix scoped, no STAR
two-pass / chimera work on this branch). Pushed to origin, up to date as of
commit `9036153`.

## Problem

User's real run: `telocal_locind` (auto-building the TElocal `.locInd` index
from `GRCm38_Ensembl_rmsk_TE.gtf`, 3.7M TE instances) had been running for
21+ hours.

## Root cause

Not the AVL/BinaryTree insertion logic (`TElocal_Toolkit.TEindex`'s own
tree code) — profiled and confirmed that's <10% of runtime. The real cost is
inside `TEfeatures.build()`'s per-line loop, in third-party
`TElocal_Toolkit` (pinned via `TElocal==1.1.3`):

1. `if ele_name not in self._elements` — membership test against a plain
   Python **list**, O(n) per check → O(n²) overall as the index grows. This
   alone measured ~98s of a ~125s run at 200k instances; at 3.7M instances
   it's the dominant cost by far.
2. Manual character-by-character attribute parsing (minor by comparison,
   but replaced too).

## Fix

`workflow/scripts/build_telocal_index.py`: added `_fast_build()`, which
reuses `TElocal_Toolkit`'s own `Node`/`BinaryTree`/`insert()` code
**verbatim, unmodified** (zero risk to query correctness — `findOvpTE` etc.
are exactly what TElocal's real quantification step calls at count time).
Only two things changed from the original `build()`:

- attribute parsing: single regex `re.compile(r'(\w+) "([^"]*)"')` instead
  of the manual loop.
- membership test: an auxiliary `set()` for O(1) lookups, while
  `te_idx._elements` is still appended to as a list in the same original
  order (required for pickle-compatibility with everything downstream that
  assumes it's a list).

## Verification (real data, not just synthetic)

- Synthetic (200k instances): identical `_nameIDmap`/`_elements`/`_length`,
  0 mismatches across 20,000 `findOvpTE` queries vs. the original `build()`.
- Real mouse data (3.7M instances, the user's actual
  `GRCm38_Ensembl_rmsk_TE.gtf`, 476MB): built via both the real pinned
  `TElocal==1.1.3` package's original `build()` and `_fast_build()`,
  identical pickled output, 0 mismatches across 50,000+ spatial queries.
  ~21h → ~1 minute.
- Cross-checked against a real pre-built index the user had downloaded
  (`GRCm38_GENCODE_rmsk_TE.gtf.locInd.gz`) — initially looked like 5,017/
  50,000 mismatches, traced to a chromosome-naming convention difference
  between the two source GTFs (`chr1`/`chrX` vs bare `1`/`X`), not a bug;
  after normalizing chromosome names and comparing by TE name, 0/50,000
  mismatches.
- Isolated venv used for all real-package testing (main conda env doesn't
  have `TElocal_Toolkit` installed — pre-existing, unrelated gap):
  `python3 -m venv` + `pip install TElocal==1.1.3 pyyaml`, since the main
  conda envs refuse global installs (PEP 668).

## Config toggle added

`telocal.indexer: fast | legacy` (default `fast`), one key inside the
existing `telocal:` block — no new config file, matching the precedent set
by `star.two_pass` (which lives on a different branch,
`feature/chimera-assembly` — **not present on this branch**, see below).

- `fast` (default): `_fast_build()`, as verified above.
- `legacy`: calls `TElocal_Toolkit`'s own unmodified `TEfeatures.build()` —
  escape hatch, since `_fast_build()` depends on that library's internal
  attributes (`_elements`, `_nameIDmap`, `_length`, `indexlist`,
  `TEindex_BINSIZE`) rather than its public API. Re-verify before bumping
  the `TElocal` version pin in `workflow/environment.yaml` /
  `workflow/default-config/versions.yaml`.

Wired through:
- `workflow/scripts/build_telocal_index.py` — new `--indexer {fast,legacy}`
  CLI flag.
- `workflow/rules/telocal.smk` — `telocal_locind` rule reads
  `config["telocal"].get("indexer", "fast")`.
- `workflow/schemas/config.schema.yaml`, `workflow/default-config/telocal.yaml`,
  `config/config.example.yaml` — schema + defaults + example, all with the
  same explanation.
- `.github/workflows/ci.yml` — guards 22/23 (dry-run: default renders
  `--indexer fast`, override renders `--indexer legacy`) and 23/23 (real
  execution: builds the index both ways against the test fixture,
  compares pickled `_nameIDmap`/`_elements`/`_length` for exact equality).
- `workflow/default-config/resources.yaml` — `telocal_locind` runtime
  budget cut from 1440 (24h) to 180 (3h); comment explains the real root
  cause and that this isn't yet re-benchmarked at full production scale,
  so there's headroom.

## telocal_locations — TE key → genomic coordinates

Separate follow-up, prompted by the user inspecting real TElocal output on
their HPC run
(`~/projects/20260805_ZPF68/.../results/telocal/2C_WT_09.cntTable.gz`):
gene rows appeared first (`head`), TE-locus rows only showed up at the
bottom (`tail`) — normal (TElocal's own output, gene rows sorted before TE
rows). But the TE rows (e.g. `hAT-N1_Mam_dup95:hAT-N1_Mam:hAT:DNA`) carry no
genomic coordinates — also normal, that's simply how TElocal's own cntTable
is shaped; it was never a pipeline bug.

Added `workflow/scripts/telocal_locations.py` + rule `telocal_locations` in
`workflow/rules/telocal.smk`: one streaming pass over `TE_GTF` (no
index/tree), emitting `results/telocal/telocal_locations.bed` — every TE
key (`transcript_id:gene_id:family_id:class_id`) mapped to its
`chrom`/`start`/`end`/`strand`, so cntTable rows can be joined against real
coordinates. Always built when `telocal.enabled: true` (like
`telocal_summary`), independent of whether `locind` is auto-built or
user-supplied.

**Follow-up: converted to BED6.** Originally a custom TSV
(`TE\tchrom\tstart\tend\tstrand`, key first, header row,
`.tsv.gz`). User asked whether this should have been BED format instead —
correct call: joining on a named column (BED's `name`, 4th column) costs
nothing over joining on column 1 (`pandas.merge(left_on=..., right_on=...)`
and `join -1 1 -2 4` are both the same plain O(n) join), so BED6 is a
strict improvement — same join cost, plus native `bedtools`/IGV
compatibility. Converted to `chrom, start0, end, name=TE key, score=".",
strand` (0-based half-open), matching the exact conventions the repo's own
`annotation_to_bed.py` already uses for `genes.bed`/`exons.bed`/`te.bed` —
no header line, `"."` score placeholder. Renamed to
`results/telocal/telocal_locations.bed`, **uncompressed** — user caught
that gzip would have silently broken the IGV-loadable half of the reason
to switch to BED at all: IGV needs plain text or `bgzip`+`tabix`, not
ordinary gzip (unlike the sibling `.bed` files' rationale, which are small
enough that compression was never worth it anyway; this file can be
millions of rows for a real genome, so the size/compatibility tradeoff was
worth being explicit about in the README).

A full-pipeline survey (Explore agent) confirmed no other output is a good
BED candidate: junction/chimera tables have 2+ distinct point-breakpoints
per row plus many essential annotation columns, count matrices have no
per-row coordinates at all, and `merged_SJ.out.tab` is contractually
STAR's own native format (fed back into STAR for cohort two-pass mapping)
and must not be touched. `te.bed` (existing, from `annotation_to_bed.py`)
looks similar but is a different artifact — one row per TE *subfamily*
(bounding-box merged across all its instances), used for chimera
junction/exon overlap tests, not a substitute for `telocal_locations`'
per-*instance* coordinates.

Verification: read TElocal's own `bin/TElocal` source directly — cntTable's
TE row keys come from the index's `getElements()`/`_elements`, built as
`transcript_id:gene_id:family_id:class_id`. Built a real `.locInd` via the
isolated venv (pinned `TElocal==1.1.3`) and confirmed `getElements()`
matches `telocal_locations.py`'s output character-for-character on the same
GTF (`TE1_dup1:TE1:TE1fam:TE1class`, `.tests/resources/te_annotation.gtf`
fixture). Dry-run + real execution against `config/test.yaml` both pass.
No new CI guard needed — it's a plain rule the existing `test` job already
exercises end-to-end (unlike `telocal.indexer`, which needed explicit
dual-path guards because it's a config-level behavior switch).

Also fixed a doc gap while in there: README.md never mentioned
`telocal.indexer: fast|legacy` even though the prior commit added it to the
pipeline — now documented, plus a new "Locus coordinates" section and
output-tree entry for `telocal_locations.tsv.gz`. Flowchart regenerated
(`telocal_locations` has no rule dependencies either side, so it renders in
its own "Other" subgraph).

## Status at end of this session

- Commit `9036153` pushed to `origin/fix/fast-telocal-indexer` (2 commits
  this session: `ce6f221` telocal.indexer toggle, `9036153`
  telocal_locations + README/agent-log updates).
- No PR opened yet for this branch (unlike `feature/chimera-assembly`,
  which has PR #1) — consider opening one for real CI validation.
- User's real project config
  (`/Volumes/jqd731/tmp/config.yaml` locally / on HPC at their `input/`
  path) was updated with `telocal.indexer: fast` under the `telocal:`
  block — that file is gitignored (`input/config.yaml`), so it's not part
  of this repo; the change lives only on their machine(s).
- `CLAUDE.md` in the repo root is intentionally **untracked** (user
  preference — do not `git add` it).

## To resume on another machine

```bash
git fetch origin
git checkout fix/fast-telocal-indexer
git pull
```

Env doesn't need changes — same pinned `TElocal` version, just the new
script/rule/config. If pointing at a fresh clone's `input/config.yaml`,
make sure `telocal.indexer: fast` (or omit it — `fast` is now the default)
is set if you want the speedup; `legacy` is the fallback.

Open threads / possible next steps:
- Decide whether to open a PR for `fix/fast-telocal-indexer` (into `main`
  or into `feature/telocal`?) to get real CI (Linux, full env) validation.
- STAR two-pass alignment (`star.two_pass: none|per_sample|cohort`) and the
  chimera-assembly work live on `feature/chimera-assembly`, not here —
  merge/rebase planning between these branches is still open.
- Not yet re-benchmarked the fast indexer at full HPC scale end-to-end
  inside the actual Snakemake run (only verified via the standalone script
  against the real GTF) — worth confirming the full `telocal_locind` rule
  behaves the same when run for real on HPC.
