#!/usr/bin/env python3
"""Deterministically generates the tiny synthetic dataset under .tests/,
used by config/test.yaml for a fast end-to-end smoke test of the workflow
(mirrors nf-core's `-profile test` pattern).

Not part of the workflow DAG -- run manually if you ever need to
regenerate/resize the test data:

    python3 .tests/generate_test_data.py

Layout produced:
    .tests/resources/genome.fa           one 50kb contig "chrT"
    .tests/resources/genome.gtf          one gene, "chrT" 10001-30000 (+)
    .tests/resources/te_annotation.gtf   one TE,   "chrT" 35001-45000 (+)
    .tests/reads/{sample}_R{1,2}.fastq.gz
    .tests/samples.csv

Genome size note: earlier versions of this script used a 10kb genome, which
turned out to sit in STAR's known small-genome crash zone (custom
references under ~10-20kb can trigger "double free or corruption" while
writing the SAindex, even with --genomeSAindexNbases correctly downscaled --
see e.g. alexdobin/STAR issues #2158, #965, #1472, #350). 50kb plus an
explicit --genomeChrBinNbits override (config/test.yaml) avoids that.
"""
import gzip
import math
import random
from pathlib import Path

HERE = Path(__file__).parent
GENOME_LEN = 50_000
CONTIG = "chrT"
GENE_START, GENE_END = 10_001, 30_000  # 1-based, inclusive (GTF convention)
TE_START, TE_END = 35_001, 45_000
READ_LEN = 50
FRAGMENT_LEN = 200  # insert size for paired-end reads

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def make_genome(seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))


def write_fasta(path, name, seq, width=70):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width] + "\n")


def write_gtf(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(row) + "\n")


def sample_read_pairs(genome, region_start, region_end, n_pairs, seed):
    """Simulate n_pairs FR paired-end reads (exact substrings, no errors)
    from within [region_start, region_end) of the + strand of `genome`
    (0-based half-open)."""
    rng = random.Random(seed)
    pairs = []
    max_start = region_end - FRAGMENT_LEN
    if max_start <= region_start:
        raise ValueError("region too small for the requested fragment length")
    for _ in range(n_pairs):
        frag_start = rng.randint(region_start, max_start)
        frag = genome[frag_start : frag_start + FRAGMENT_LEN]
        r1 = frag[:READ_LEN]
        r2 = revcomp(frag[-READ_LEN:])
        pairs.append((r1, r2))
    return pairs


def write_fastq_gz(path, reads, name_prefix):
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(reads):
            qual = "I" * len(seq)
            fh.write(f"@{name_prefix}_{i}/1\n{seq}\n+\n{qual}\n")


def main():
    resources_dir = HERE / "resources"
    reads_dir = HERE / "reads"
    resources_dir.mkdir(parents=True, exist_ok=True)
    reads_dir.mkdir(parents=True, exist_ok=True)

    genome = make_genome(seed=1)
    write_fasta(resources_dir / "genome.fa", CONTIG, genome)

    write_gtf(
        resources_dir / "genome.gtf",
        [
            (
                CONTIG,
                "test",
                "exon",
                str(GENE_START),
                str(GENE_END),
                ".",
                "+",
                ".",
                'gene_id "GENE1"; transcript_id "GENE1"; gene_name "GENE1";',
            )
        ],
    )

    # TEtranscripts' curated TE GTF format: "exon" feature, gene_id/
    # transcript_id/family_id/class_id attributes (matches the format of
    # the authors' official prebuilt TE GTFs, e.g. mm10/hg38 rmsk.gtf).
    write_gtf(
        resources_dir / "te_annotation.gtf",
        [
            (
                CONTIG,
                "test",
                "exon",
                str(TE_START),
                str(TE_END),
                ".",
                "+",
                ".",
                'gene_id "TE1"; transcript_id "TE1_dup1"; '
                'family_id "TE1fam"; class_id "TE1class";',
            )
        ],
    )

    # 4 samples, 2 conditions x 2 reps, all paired-end. Treatment samples
    # get more gene-region reads and fewer TE-region reads than control
    # (and vice versa) purely so DESeq2 sees non-zero variance between
    # groups -- there's no biological meaning to the direction/magnitude,
    # this is a connectivity smoke test, not a validation dataset.
    sample_plan = {
        "treatment_rep1": {"condition": "treatment", "gene_pairs": 160, "te_pairs": 40, "seed": 101},
        "treatment_rep2": {"condition": "treatment", "gene_pairs": 150, "te_pairs": 45, "seed": 102},
        "control_rep1": {"condition": "control", "gene_pairs": 80, "te_pairs": 100, "seed": 201},
        "control_rep2": {"condition": "control", "gene_pairs": 90, "te_pairs": 95, "seed": 202},
    }

    # 0-based half-open coordinates for read sampling (GTF above is 1-based
    # inclusive; convert once here).
    gene_region = (GENE_START - 1, GENE_END)
    te_region = (TE_START - 1, TE_END)

    rows = ["sample,fastq_1,fastq_2,strandedness,condition"]
    for sample, plan in sample_plan.items():
        gene_pairs = sample_read_pairs(
            genome, *gene_region, plan["gene_pairs"], seed=plan["seed"]
        )
        te_pairs = sample_read_pairs(
            genome, *te_region, plan["te_pairs"], seed=plan["seed"] + 1
        )
        all_pairs = gene_pairs + te_pairs
        r1_reads = [p[0] for p in all_pairs]
        r2_reads = [p[1] for p in all_pairs]

        fq1 = reads_dir / f"{sample}_R1.fastq.gz"
        fq2 = reads_dir / f"{sample}_R2.fastq.gz"
        write_fastq_gz(fq1, r1_reads, sample)
        write_fastq_gz(fq2, r2_reads, sample)

        rows.append(
            f"{sample},.tests/reads/{sample}_R1.fastq.gz,"
            f".tests/reads/{sample}_R2.fastq.gz,auto,{plan['condition']}"
        )

    (HERE / "samples.csv").write_text("\n".join(rows) + "\n")
    print("Wrote test data to", HERE)

    # STAR's own formulas (see the STAR manual, section 2.2.5) -- used to
    # compute config/test.yaml's star.index_extra. Re-run this script and
    # copy the printed values there if you change GENOME_LEN.
    sa_index_nbases = min(14, int(math.log2(GENOME_LEN) / 2 - 1))
    chr_bin_nbits = min(18, int(math.log2(GENOME_LEN)))
    print(
        f"Recommended star.index_extra for GENOME_LEN={GENOME_LEN}: "
        f'"--genomeSAindexNbases {sa_index_nbases} '
        f'--genomeChrBinNbits {chr_bin_nbits}"'
    )


if __name__ == "__main__":
    main()
