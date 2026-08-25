# TElocal fast indexer — session log

Branch: `fix/fast-telocal-indexer` (based on `feature/telocal`, deliberately
*not* on `feature/chimera-assembly` — keeps this fix scoped, no STAR
two-pass / chimera work on this branch). Pushed to origin, up to date as of
commit `ce6f221`.

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

## Status at end of this session

- Commit `ce6f221` pushed to `origin/fix/fast-telocal-indexer`.
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
