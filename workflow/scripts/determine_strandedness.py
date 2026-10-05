"""Turn RSeQC infer_experiment.py output into a TEtranscripts --stranded value.

RSeQC infer_experiment.py reports, for single-end data:
    Fraction of reads explained by "++,--":   <forward-type fraction>
    Fraction of reads explained by "+-,-+":   <reverse-type fraction>
and for paired-end data:
    Fraction of reads explained by "1++,1--,2+-,2-+": <forward-type fraction>
    Fraction of reads explained by "1+-,1-+,2++,2--": <reverse-type fraction>

Mapping to TEtranscripts/TEcount's --stranded {no,forward,reverse}:
  forward-type fraction >= min_fraction -> "forward" (e.g. QIAseq stranded / "second-strand")
  reverse-type fraction >= min_fraction -> "reverse" (e.g. Illumina TruSeq stranded / "first-strand")
  neither clears min_fraction          -> "no" (the operative --stranded value; see below)

Writes TWO files:
  output.txt   the operative --stranded value (no/forward/reverse) -- what
               TEcount/TElocal actually run with for this sample (when its
               effective mode is "auto").
  output.call  a richer, report-only label (forward/reverse/no/undetermined):
               below min_fraction, the dominant fraction is checked against
               a second, lower threshold (balanced_max) -- under it, the
               library is confidently unstranded ("no"); between the two,
               neither confidently stranded nor confidently balanced, so
               it's reported as "undetermined" rather than silently folded
               into "no". Quantification still uses "no" for undetermined
               samples (the safe choice -- it discards no reads), but
               strandedness_check_mqc.py surfaces the distinction.
"""

import re
import sys

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input.txt) as fh:
    text = fh.read()

fractions = {}
for line in text.splitlines():
    m = re.search(r'explained by "([^"]+)":\s*([\d.]+)', line)
    if m:
        pattern, value = m.group(1), float(m.group(2))
        fractions[pattern] = value

if not fractions:
    # RSeQC may produce no parseable output for an empty BAM (zero reads).
    # Default to "no" (unstranded) rather than crashing the pipeline --
    # the downstream tools will simply find zero counts, which is correct
    # for a sample with no data.
    print(
        f"WARNING: could not parse any 'explained by' fractions from "
        f"{snakemake.input.txt} (empty BAM?). Defaulting to 'no' (unstranded).",
        file=sys.stderr,
    )
    with open(snakemake.output.txt, "w") as fh:
        fh.write("no\n")
    with open(snakemake.output.call, "w") as fh:
        fh.write("no\n")
    sys.exit(0)

forward_value = 0.0
reverse_value = 0.0
for pattern, value in fractions.items():
    first_token = pattern.split(",")[0]  # "++" / "1++" / "+-" / "1+-"
    if first_token.endswith("++"):
        forward_value = value
    elif first_token.endswith("+-"):
        reverse_value = value

min_fraction = float(snakemake.params.min_fraction)
balanced_max = float(snakemake.params.balanced_max)

print(f"Parsed fractions: {fractions}", file=sys.stderr)
print(
    f"forward-type fraction: {forward_value:.4f}, "
    f"reverse-type fraction: {reverse_value:.4f}, "
    f"min_fraction threshold: {min_fraction}, balanced_max threshold: {balanced_max}",
    file=sys.stderr,
)

if forward_value >= min_fraction:
    call = "forward"
elif reverse_value >= min_fraction:
    call = "reverse"
elif max(forward_value, reverse_value) < balanced_max:
    call = "no"  # confidently unstranded
else:
    call = "undetermined"  # skewed, but not enough to call confidently

# The operative --stranded value TEcount/TElocal run with: undetermined
# collapses to "no", the safe choice (it discards no reads, only loses
# directionality signal).
stranded = "no" if call == "undetermined" else call

print(
    f"--> report call: {call}; operative --stranded value: {stranded}",
    file=sys.stderr,
)

with open(snakemake.output.txt, "w") as fh:
    fh.write(stranded + "\n")
with open(snakemake.output.call, "w") as fh:
    fh.write(call + "\n")
