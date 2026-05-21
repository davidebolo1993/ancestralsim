#!/usr/bin/env python3
"""Add extra msprime-modeled divergence to selected haplotype FASTAs.

The input FASTAs are expected to contain one sequence each. Original headers are
preserved in the optional PanSN output, while the per-haplotype outputs can use a
shared gargammel-compatible header.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import msprime
import tskit


BASES = ("A", "C", "G", "T")


def read_single_fasta(path: Path) -> tuple[str, str]:
    header = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"{path} contains more than one FASTA record")
                header = line[1:].split()[0]
            else:
                seq_parts.append(line)
    if header is None:
        raise ValueError(f"{path} does not contain a FASTA header")
    return header, "".join(seq_parts).upper()


def write_fasta(path: Path, header: str, sequence: str, width: int = 80) -> None:
    with path.open("w") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            handle.write(sequence[i : i + width] + "\n")


def is_pansn(header: str, delimiter: str) -> bool:
    parts = header.split(delimiter, 2)
    return len(parts) == 3 and all(parts) and parts[1].isdigit()


def branch_tree_sequence(sequence_length: int, generations: float) -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length=sequence_length)
    parent = tables.nodes.add_row(time=generations)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(left=0, right=sequence_length, parent=parent, child=child)
    tables.sort()
    return tables.tree_sequence()


def transition_base(current: str, rng: random.Random) -> str:
    choices = [base for base in BASES if base != current]
    return rng.choice(choices)


def mutation_positions(
    sequence_length: int,
    rate: float,
    generations: float,
    seed: int,
    model: str,
) -> dict[int, str]:
    ts = branch_tree_sequence(sequence_length, generations)
    if model == "jc69":
        mutation_model = msprime.JC69()
    else:
        raise ValueError(f"Unsupported mutation model: {model}")

    mts = msprime.sim_mutations(
        ts,
        rate=rate,
        model=mutation_model,
        random_seed=seed,
        discrete_genome=True,
    )
    out: dict[int, str] = {}
    sample_id = mts.samples()[0]
    for variant in mts.variants(samples=[sample_id]):
        pos = int(variant.site.position)
        allele = variant.alleles[int(variant.genotypes[0])]
        if allele in BASES:
            out[pos] = allele
    return out


def diverge_sequence(
    sequence: str,
    rate: float,
    generations: float,
    seed: int,
    model: str,
) -> tuple[str, int, int]:
    rng = random.Random(seed + 104729)
    seq = list(sequence)
    candidate_positions = mutation_positions(len(seq), rate, generations, seed, model)
    applied = 0
    skipped_non_acgt = 0

    for pos, msprime_base in candidate_positions.items():
        current = seq[pos].upper()
        if current not in BASES:
            skipped_non_acgt += 1
            continue
        if msprime_base == current:
            msprime_base = transition_base(current, rng)
        seq[pos] = msprime_base
        applied += 1

    return "".join(seq), applied, skipped_non_acgt


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Overlay msprime-modeled divergence onto selected haplotypes."
    )
    parser.add_argument("--hap1", required=True, type=Path)
    parser.add_argument("--hap2", required=True, type=Path)
    parser.add_argument("--out-hap1", required=True, type=Path)
    parser.add_argument("--out-hap2", required=True, type=Path)
    parser.add_argument("--pansn-output", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--header", default=None)
    parser.add_argument("--rate", required=True, type=float)
    parser.add_argument("--generations", required=True, type=float)
    parser.add_argument("--model", default="jc69", choices=["jc69"])
    parser.add_argument("--seed", type=int)
    parser.add_argument("--pansn-delimiter", default="#")
    parser.add_argument("--require-pansn", action="store_true")
    args = parser.parse_args()

    if args.rate < 0:
        raise ValueError("--rate must be >= 0")
    if args.generations < 0:
        raise ValueError("--generations must be >= 0")

    header1, seq1 = read_single_fasta(args.hap1)
    header2, seq2 = read_single_fasta(args.hap2)
    if args.require_pansn:
        bad = [
            header
            for header in (header1, header2)
            if not is_pansn(header, args.pansn_delimiter)
        ]
        if bad:
            raise ValueError(
                "Input FASTA headers do not follow PanSN "
                f"sample#haplotype#contig format: {', '.join(bad)}"
            )

    seed = args.seed if args.seed is not None else random.SystemRandom().randint(1, 2**32 - 1)
    seq1_div, hap1_mutations, hap1_skipped = diverge_sequence(
        seq1, args.rate, args.generations, seed, args.model
    )
    seq2_div, hap2_mutations, hap2_skipped = diverge_sequence(
        seq2, args.rate, args.generations, ((seed + 1000003 - 1) % (2**32 - 1)) + 1, args.model
    )

    output_header1 = args.header if args.header else header1
    output_header2 = args.header if args.header else header2
    write_fasta(args.out_hap1, output_header1, seq1_div)
    write_fasta(args.out_hap2, output_header2, seq2_div)

    if args.pansn_output is not None:
        with args.pansn_output.open("w") as handle:
            for header, sequence in ((header1, seq1_div), (header2, seq2_div)):
                handle.write(f">{header}\n")
                for i in range(0, len(sequence), 80):
                    handle.write(sequence[i : i + 80] + "\n")

    if args.report is not None:
        callable1 = sum(base in BASES for base in seq1)
        callable2 = sum(base in BASES for base in seq2)
        expected1 = args.rate * args.generations * callable1
        expected2 = args.rate * args.generations * callable2
        with args.report.open("w") as handle:
            handle.write("haplotype\tlength\tcallable_acgt\texpected_mutations\tapplied_mutations\tskipped_non_acgt\tseed\n")
            handle.write(
                f"{header1}\t{len(seq1)}\t{callable1}\t{expected1:.6f}\t"
                f"{hap1_mutations}\t{hap1_skipped}\t{seed}\n"
            )
            handle.write(
                f"{header2}\t{len(seq2)}\t{callable2}\t{expected2:.6f}\t"
                f"{hap2_mutations}\t{hap2_skipped}\t{((seed + 1000003 - 1) % (2**32 - 1)) + 1}\n"
            )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
