#!/usr/bin/env python3
import re
import sys

infer_file = sys.argv[1]
output_file = sys.argv[2]
min_fraction = float(sys.argv[3])

with open(infer_file) as fh:
    text = fh.read()

fractions = {}
for line in text.splitlines():
    m = re.search(r'explained by "([^"]+)":\s*([\d.]+)', line)
    if m:
        pattern, value = m.group(1), float(m.group(2))
        fractions[pattern] = value

if not fractions:
    raise ValueError(
        f"Could not parse any 'explained by' fractions from {infer_file}."
    )

forward_value = 0.0
reverse_value = 0.0
for pattern, value in fractions.items():
    first_token = pattern.split(",")[0]
    if first_token.endswith("++"):
        forward_value = value
    elif first_token.endswith("+-"):
        reverse_value = value

if forward_value >= min_fraction:
    stranded = "forward"
elif reverse_value >= min_fraction:
    stranded = "reverse"
else:
    stranded = "no"

with open(output_file, "w") as fh:
    fh.write(stranded + "\n")
