# Chimera evidence: what has actually been measured

This pipeline does not rank or score gene-TE chimera candidates (see
`workflow/scripts/chimera_evidence.py` and `chimera_candidates_table_mqc.py`).
The report's "How to weigh this evidence" section explains what each signal
*is* and how to read it, in wording that holds for any run. This document is
the other half: the project's own measurements of those signals, each tied
to the specific cohort and restrictions it was measured under, so a reader
can judge how far -- if at all -- a finding here should inform reading their
own run.

None of the measurements below has established a weighting for any signal.
They are why some signals are already flagged "Read with care" or "Not
validated" rather than treated as settled; they are not a substitute for
checking your own `candidates.tsv.gz`.

Signals are listed in the same order as the report's guide and the
Candidates table's own columns: Screens / Found by, CR motif, SJ motif,
Assembly strand, Replicated, TElocal reads, Read depth.

## Measured so far

### Screens / Found by -- cross-screen agreement near chance

- **What was measured:** whether pairs found by more than one screen were
  found disproportionately more often than chance would predict, as part of
  a dimension-by-dimension comparison of the evidence columns.
- **Cohort:** a 4-sample mouse run.
- **Restriction:** only pairs the chimeric-reads (CR) screen found --
  assembly-only and splice-junctions-only pairs were not part of this
  measurement, so it says nothing about them. This is a biased subset of
  the full candidate catalogue, not the cohort as a whole.
- **Finding:** cr+assembly agreement on that subset came out near its chance rate.
- **What it shows:** on this one biased subset, two screens agreeing on a
  pair was not the strong, near-independent confirmation the screens'
  different blind spots would suggest.
- **What it doesn't show:** anything about the SJ screen's agreement with
  either of the other two (it had not been run on real data at the time of
  this measurement), three-screen agreement, or whether the same holds on
  the full, unrestricted catalogue, on a different cohort, or on a different
  organism.
- **Why this matters here:** this result is what directly contradicted the
  top rung of the four-tier confidence ladder this pipeline used to render
  (removed; see "Why the ladder and the per-screen top-N tables were
  removed" below).

### TElocal reads -- expression looked like the opposite of support

- **What was measured:** whether the TE locus being called "expressed" by
  TElocal (`telocal_expressed`) correlated with the CR screen's own splice
  motif rate, and how common TElocal-expressed loci were among CR-screen
  pairs.
- **Cohort:** the same 4-sample mouse run.
- **Restriction:** the same biased subset -- only pairs the chimeric-reads
  screen found.
- **Finding:** 91% of those pairs had an expressed TE locus, so on its own
  the flag discriminated almost nothing within that subset. Where the locus
  was expressed, the canonical (splice-motif) rate was *lower*, not higher:
  6.7% vs 10.2% canonical, n = 19,503 events.
- **What it shows:** on this subset, an expressed TE locus did not look like
  independent support for a chimera call -- if anything, the direction was
  the opposite of what "more expression, more support" would predict.
- **What it doesn't show:** causation, or whether the same direction holds
  outside this one biased subset. It is mechanistically unsurprising either
  way -- a highly expressed locus yields more reads and so both more chances
  for a genuine chimera to be observed and more chances for a template-
  switching artifact -- so this result alone is not enough to demote the
  signal.
- **Related note:** `chimera_evidence.py`'s own docstring additionally
  describes a threshold-based version of this same comparison
  (`telocal_count > 10` vs `<= 10`) rather than the expressed/not-expressed
  split above; both come from the same run and the same general finding
  (lower canonical rate at higher TElocal signal), but were computed as two
  different splits of the same underlying data and have not been
  reconciled into one number.

## Not yet measured

- **SJ-screen agreement** with either of the other two screens (CR motif's
  chance-rate comparison above predates the SJ screen running on real data).
- **Three-screen agreement** (`n_screens == 3`).
- **Condition-aware replication** -- whether `Replicated` (seen in more than
  one sample) means something different when the samples span different
  experimental conditions vs. technical replicates of the same one.

## Why the ladder and the per-screen top-N tables were removed

This pipeline used to rank chimera candidates two different ways, both
removed:

- A **four-tier confidence ladder**, the report's first chimera section.
  Removed because no experiment in this project established the relative
  weight of its rungs, and the Screens / Found by finding above directly
  contradicted its top rung. The same class of result had already removed
  TElocal expression from the ladder before that (commit `0d04e43`).
- **Three competing per-screen top-N tables** (one per screen, each with its
  own ranking key), removed in commit `d927c8f` for the same reason: no
  validated weighting, and three different orderings in one report read as
  three different verdicts.

Ranking on an unvalidated weighting is worse than not ranking, because a
rank or tier in a report is read as a verdict. The Candidates table
(`chimera_candidates_table_mqc.py`) is the replacement: every pair, every
evidence column, sorted by Screens only and re-sortable by clicking any
other header -- a real count the reader can act on, not a score the
pipeline asserts.

## Next step

None of the above is enough to change how any signal is weighted, used as a
filter, or combined with any other. The next step is a dedicated analysis on
real project data -- ideally without the chimeric-reads-screen restriction
that limited every measurement above -- before any of it goes into the
pipeline's own logic rather than this document.
