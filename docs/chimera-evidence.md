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

### TE orientation -- LTR-initiated calls are strongly sense-biased

- **What was measured:** the TE's strand relative to its gene (sense /
  antisense), per screen x chimera type x TE class, against a background
  of TE copies of the same class at the same position relative to the
  same genes (upstream / downstream in distance rings, intronic, exonic),
  combined across positions with a Mantel-Haenszel odds ratio. First on a
  4-sample run, then repeated on the full 84-sample cohort with the
  current labels (strand rule, `chimera_status`, assembly TE choice).
- **Cohort:** a mouse oocyte/embryo series (unstranded, single-end): 4
  samples, then all 84 (WT and mKO at GV, zygote, 2-cell, 8-cell,
  blastocyst).
- **Finding (84 samples):** LTR `te_initiated` is strongly sense-biased --
  SJ 66.5% vs 43.5% expected (OR 2.6, n = 12,620), assembly 71.1% vs
  44.2% (OR 3.2, n = 7,608) -- and the bias rises with screen agreement
  (SJ 60% -> 79%, assembly 63% -> 79% for 1 vs >= 2 screens) and with
  replication (SJ 57% -> 72%, assembly 63% -> 83%). On the 4-sample run
  the same held (SJ 69.9% vs 43.7%, assembly 66.8% vs 43.0%). Chimeric-
  reads LTR `te_initiated` sits near background (53.1% vs 49.1%, OR 1.17,
  n = 565), but only for pairs that screen alone finds (46% vs 49%); the
  ones another screen also finds are 61% vs 49% sense. LINE `te_initiated`
  leans antisense after position matching (SJ 33.4% vs 40.3%, OR 0.74;
  assembly 32.9% vs 42.4%, OR 0.66) -- on 4 samples this looked mostly
  like a position effect. SINE `te_exonized` is antisense-biased (SJ 27.9%
  vs 45.8%, OR 0.46, n = 11,055), the mouse counterpart of antisense Alu
  exonization. LTR `te_exonized` is only weakly sense-biased (OR 1.1-1.2);
  the intronic-LTR-promoter signal sits in `te_initiated` (internal) and
  `annotated_promoter_embedded_te` (LTR OR 2.7) instead. LTRs in
  `antisense_to_gene` calls are ~75% sense to the antisense transcript
  itself (22-26% sense to the gene vs ~47% expected) -- the pattern of
  LTR-driven antisense transcription, not of background.
- **What it shows:** orientation separates real LTR-driven initiation from
  background in the SJ and assembly screens, and screen agreement and
  replication pick out more of it. A chimeric-reads call on its own looks
  like background; confirmed by another screen, it does not.
- **What it doesn't show:** a validated weighting, or which individual
  calls are real. One cohort, pairs not independent -- p-values are
  descriptive.

### Classification problems found on the same run (fixed)

Measured with a one-off check of every call's geometry and annotation:
- ~500 SJ and ~400 assembly `te_initiated` calls had their TE downstream of
  the whole gene (and ~400 assembly `te_terminated` calls had it upstream):
  antisense transcription through the gene's exon. Now typed
  `antisense_to_gene`; after the fix no wrong-side calls remained.
- 25.1% of SJ gene-TE junctions were annotated GTF introns (STAR's own
  `annotated` flag said 51.1% -- unreliable under `star.two_pass: cohort`).
  Now typed `annotated_splice`.
- 42.8% of assembly `te_terminated` calls had their TE in the gene's
  annotated last exon (an ordinary 3' UTR TE). Now typed
  `annotated_terminal_exon_embedded_te`.
- 97.5% of chimeric-read gene-TE events were trans or >200 kb from the
  gene. They no longer count toward candidates
  (`chimera.chimeric_reads.max_gene_te_distance`).

### Replicated / SJ samples -- cohort 2-pass inflates them slightly

- **What was measured:** the same 4 samples run twice at one pipeline
  commit, differing only in `star.two_pass` (`cohort`: junctions pooled from
  every sample's first pass are inserted into one index; `per_sample`: each
  sample's own). Compared: `sj_max_samples` for SJ chimera calls, the
  `multi_sample` (Replicated) rate, candidate overlap and the Screens = 3 set.
- **Cohort:** the same 4-sample mouse oocyte run (unstranded, single-end).
- **Restriction:** one cohort of 4 samples, all one stage; with more samples
  the pooled index holds more junctions, so the effect may grow.
- **Finding:** of 11,664 pairs both runs called by SJ, 864 had a higher
  `sj_max_samples` under cohort mode and 22 a lower one; the share seen in
  all 4 samples was 2,246 vs 1,814 pairs. Replicated: 43.3% vs 42.2% of
  SJ-found pairs (32.9% vs 31.8% overall); 121 shared pairs were replicated
  only under cohort mode, 18 only under per-sample. Candidates otherwise
  agreed closely: 15,318 shared, 146 cohort-only, 329 per-sample-only
  (Jaccard 0.97); Screens = 3 was 21 vs 22 pairs, all 21 shared. STAR's own
  `annotated` flag marked 79-82% of SJ.out.tab junctions under cohort mode
  and 99.7-99.8% under per-sample mode, so it is unusable as "annotated in
  the GTF" in either mode.
- **What it shows:** cohort 2-pass lets a junction seen in one sample's first
  pass be detected in others on the second pass, which lifts sample counts:
  a pair's `sj_max_samples` and Replicated can partly reflect the pooled
  index rather than independent detection. The size here is small (about
  one percentage point of Replicated) but one-directional.
- **What it doesn't show:** whether the extra detections are real junctions
  that per-sample mode misses or index-induced artifacts, or how the effect
  scales with cohort size.
- **Why this matters here:** Replicated and SJ samples are read as
  recurrence across independent libraries; under `star.two_pass: cohort`
  they are not fully independent.

### Full cohort -- what the classification fixes changed (v0.14.0 vs pre-release 0.15.0, commit 80f1ac6)

- **What was measured:** the same 84 samples run end to end with v0.14.0
  (before direction-based typing, the strand rule, the known-structure
  classes, the far/trans chimeric-read filter and the SJ screen) and with
  a pre-release 0.15.0 build, commit 80f1ac6 (SJ screen on, cohort
  2-pass); the released 0.15.0 adds the changes recorded in the next
  sections. Compared: candidate pairs, where
  the old pairs went, the top-Screens set, chimera types per screen, and
  the TEcount / TElocal count matrices.
- **Cohort:** a mouse oocyte/embryo knockout series, 84 samples (WT
  and mKO at GV, zygote, 2-cell, 8-cell, blastocyst), unstranded single-end.
- **Restriction:** one cohort; no per-sample 2-pass run, so cohort 2-pass
  effects on Replicated cannot be separated from cohort size here.
- **Finding:** candidates fell from 454,592 to 79,915 pairs. Of the 423,943
  old pairs no longer present, 405,990 (96%) were chimeric-read pairs on
  another chromosome or >200 kb from the gene, 11,817 (almost all
  assembly-only) are now typed only as antisense or known gene structure,
  and 6,136 are no longer seen. The SJ screen supplies most new candidates
  (62% are SJ-only); the chimeric-reads screen contributes to 4,620 pairs.
  310 pairs are found by all three screens; 1,310 of the 1,312 pairs v0.14.0
  found with both of its screens are still candidates. Quantification barely
  moved: per-sample Spearman old vs new >= 0.994 (TEcount) and >= 0.954
  (TElocal), totals within 0.03%, TE counts ~0.9% lower (chimeric
  supplementary records are no longer counted twice). Replicated was 52.9%
  of SJ-found pairs. STAR's own `annotated` flag marked 46.1% of typed SJ
  junctions, none of which is a GTF intron.
- **What it shows:** on a full cohort, most of the old candidate catalogue
  was random-partner chimeric-read noise and known gene structure; the
  remaining candidates rest mostly on the SJ and assembly screens.
- **What it doesn't show:** which of the remaining candidates are real, or
  how the numbers would change with per-sample 2-pass or on another cohort.
- **Why this matters here:** results from v0.14.0 or earlier are not
  comparable pair-for-pair with 0.15.0; rerun before comparing candidates.

### Full cohort -- counting annotated and antisense chimeras

- **What was measured:** the same 84-sample run rebuilt after every
  typed gene-TE chimera started counting toward Screens, with
  `chimera_status` saying which kind (novel / annotated / antisense).
  Compared: candidate pairs, the status breakdown, the top-Screens set, and
  which screens make antisense calls.
- **Cohort:** as above (84 samples, unstranded single-end).
- **Restriction:** one cohort; unstranded, so the chimeric-reads screen
  cannot call antisense here.
- **Finding:** candidates rose from 79,915 to 101,611 pairs; none were
  lost. The 21,696 new pairs are 11,189 antisense, 10,107 annotated and 400
  both, found by the assembly (11,062), SJ (8,505) or both screens (2,129).
  Of the old pairs, 7,702 gained an annotated or antisense kind and 3,847
  gained a screen. Pairs found by all three screens rose from 310 to 803;
  of these only 170 are `novel` alone (336 `novel+annotated`, 222
  `novel+annotated+antisense`, 75 `novel+antisense`). Antisense calls came
  from SJ (8,796 pairs) and assembly (6,260) only; the chimeric-reads
  screen made none. Across all pairs: 72,213 `novel`, 11,189 `antisense`,
  10,107 `annotated`, 5,821 `novel+annotated`, 1,363 `novel+antisense`,
  518 `novel+annotated+antisense`, 400 `annotated+antisense`.
- **What it shows:** many pairs that the three screens agree on are partly
  known TE-driven transcripts or antisense transcription, and their
  per-screen read counts include those reads. Known examples now back as
  candidates: a TE-in-3'-UTR pair a DE analysis on the assembly matrix had
  ranked first (`annotated`), and an annotated TE promoter also seen as
  novel splicing (`novel+annotated`, all three screens).
- **What it doesn't show:** whether the annotated or antisense kinds
  differ in how often they are real, or how a stranded library would split
  the chimeric-reads screen's novel calls.
- **Why this matters here:** for new chimeras, filter on Status =
  `novel`; for every TE-driven transcript, use all statuses. Screens alone
  no longer separates the two.

### Assembly TSS-in-TE gate -- what it rejects looks like background

- **What was measured:** on the full cohort, the assembly classifier was
  re-run with `chimera.assembly.require_tss_in_te` off. For every
  transcript the gate rejects (a TE in the first exon, but the TSS more
  than `breakpoint_tolerance` = 5 bp outside it): the TSS-to-TE distance,
  whether the TE is sense to the transcript, and whether the
  chimeric-reads or SJ screen also found that gene-TE pair. Baselines:
  calls that pass the gate, and annotated first exons overlapping a TE.
- **Cohort:** as above (84 samples, unstranded single-end).
- **Restriction:** one cohort; StringTie merged assembly only.
- **Finding:** the gate rejects 17,378 transcripts (with it off they would
  be 13,889 `te_initiated` / `te_initiated_intergenic` and 3,489
  `annotated_promoter_embedded_te`). Calls passing the gate are 65% sense
  and 51% (`te_initiated`) / 42% (annotated promoter) found by another
  screen. Rejected transcripts are 54-60% sense at every distance (6-50,
  51-100, 101-300, 301-1000 bp) -- close to the 48-54% of annotated first
  exons with a TE -- and only 20-26% are found by another screen. The
  nearest distances (6-50 bp) have the lowest confirmation (20%); a 300 bp
  tolerance would add ~15,000 transcripts, ~80% unconfirmed.
- **What it shows:** imprecise StringTie 5' ends are not hiding many TE
  promoters behind the gate; what it rejects looks like ordinary 5' UTRs
  that contain a TE. Individual real cases exist (an MTB LTR promoter
  281 bp downstream of an assembled TSS, found by both other screens).
- **What it doesn't show:** whether a TSS-aware assembler or 5' data would
  recover more.
- **Why this matters here:** `require_tss_in_te` stays on with the
  ordinary tolerance; a TE promoter missed by the assembly is still found
  by the junction-based screens and counted in their matrices.

Measured on the same run: where several TEs hit one exon, the assembly
classifier named the one with the lowest coordinate. Of 2,911 initiation
calls with several TEs at the TSS, 1,250 named a different TE than the one
the first exon runs through, 969 of them a TE with <= 5 bp in that exon
(e.g. an L1 ending at an MT2_Mm promoter's TSS). The classifier now names
the TE at the TSS (first exon), with most bases in the exon (internal
exon), or nearest the splice acceptor (last exon). No call changes type;
9,599 of 158,434 classified transcripts name a different TE, 4,517 of them
a different TE class. Of the renamed calls, the share whose gene-TE pair
the chimeric-reads or SJ screen also found, old pick -> new pick:
`te_initiated` 4.7% -> 50.4% (381 calls), `annotated_promoter_embedded_te`
5.7% -> 46.6% (88), `te_exonized` 51.9% -> 55.0% (545), `te_terminated`
6.4% -> 41.3% (2,267), `annotated_terminal_exon_embedded_te` 18.7% ->
46.2% (4,585), `antisense_to_gene` 11.1% -> 33.8% (769). For last exons,
naming the TE at the transcript's 3' end instead was tried and rejected:
confirmation of the renamed terminal calls fell to 7-15%. What makes a
transcript TE-terminated is the splice into a TE-derived last exon, and
StringTie's 3' ends are imprecise.

### Assembly last-exon TE distance -- termination signal sits near the acceptor

- **What was measured:** for assembly `te_terminated` and
  `annotated_terminal_exon_embedded_te` calls, the distance from the last
  exon's splice acceptor to the named TE, against the TE's orientation to
  the transcript (a TE that ends a transcript must be sense to it to
  supply its polyadenylation signal). Other-screen confirmation is not
  used here: the junction screens only see a TE where a junction lands, so
  it would favour distance 0 by construction.
- **Cohort:** as above (84 samples, unstranded single-end).
- **Restriction:** one cohort; StringTie 3' ends are imprecise, so whether
  the TE also holds the transcript's 3' end was not informative (29-37% at
  every distance).
- **Finding:** `te_terminated` (n = 13,633): at 0-5 bp (the TE takes the
  splice, 39% of calls) 43% sense overall but SINEs only 28% -- antisense
  SINEs supplying an acceptor, the exonization pattern; at 6-50 bp 68%
  sense and 51-200 bp 63% (SINEs ~70%) -- 16% of calls; at 201-500 bp 50%;
  beyond 500 bp 42-46%, LTRs 27-38% -- at or below background (~43-51%),
  46% of calls. `annotated_terminal_exon_embedded_te` (n = 19,565) has the
  same shape: 69-76% sense at 6-200 bp, 51-59% beyond 200 bp.
- **What it shows:** the TE-terminated signal is concentrated within ~200
  bp of the last exon's acceptor; farther out the TE is mostly UTR content.
  The 0 bp group is real but a different mechanism (TE-derived terminal
  exons).
- **What it doesn't show:** a validated cut-off.
- **Why this matters here:** rather than a gate, the distance is reported
  as `te_acceptor_distance_bp` (assembly transcripts) and
  `assembly_te_acceptor_distance_bp` / "Assembly last-exon TE distance"
  (candidates), to sort or filter on.

### "Other" TE classes -- mostly zinc-finger coding repeats, sense by construction

- **What was measured:** the candidates whose TE is not LTR / LINE / SINE /
  DNA (RepeatMasker classes Unknown, Satellite, Other, RNA, RC and the
  "?" classes), by family, gene, gene biotype, TE position and orientation,
  and whether the TE overlaps the gene's own CDS (reference GTF).
- **Cohort:** the same 4-sample mouse oocyte run, at the pipeline commit
  that keeps only chimera calls (359 such pairs; ~1,260 before the
  classification fixes above).
- **Restriction:** one cohort; one annotation (GENCODE vM23, the
  TEtranscripts GRCm38 rmsk GTF).
- **Finding:** one family dominates: MurSatRep1 (class Unknown), 183 of the
  359 pairs, 98% sense. 162 of those 183 genes are zinc-finger genes
  (Zfp*, Zscan*, Gm* in KZFP clusters), and in 95 pairs the "TE" overlaps
  the gene's own CDS -- all 95 sense. The rest of the class is small:
  Satellite 73 pairs (74% sense), Other / RMER1 66 (45-62% sense), RNA 7.
- **What it shows:** the near-100% sense rows are RepeatMasker annotating
  the zinc-finger coding repeats themselves as a repeat family, so the TE
  lies on the gene's strand by construction. They are gene structure, not
  a TE driving the gene.
- **What it doesn't show:** whether any of these pairs is a real chimera;
  a call on such a pair is not excluded, only uninformative for
  orientation.
- **Why this matters here:** read TE orientation for the Unknown /
  MurSatRep1 rows (and any TE overlapping the gene's CDS) as annotation,
  not evidence; keep these classes out of orientation summaries, as the
  orientation measurement above did.

## Not yet measured

- **SJ-screen agreement** as such (the orientation result above shows
  agreement tracking real signal, but agreement was not measured against
  a chance rate the way the CR+assembly comparison was).
- **Three-screen agreement** (`n_screens == 3`).
- **Any of the above on a second cohort.**
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
