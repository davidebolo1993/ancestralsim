#!/usr/bin/env python3
"""Add extra msprime-modeled divergence to selected haplotype FASTAs.

Two modes are available:
  terminal: add mutations along an extra terminal branch.
  split: simulate a modern/ancient population split and overlay the resulting
         ancient-vs-modern differences onto the selected haplotypes.

Original headers are preserved in the optional PanSN output, while the
per-haplotype outputs can use a shared gargammel-compatible header.
"""

from __future__ import annotations

import argparse
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import msprime
import tskit


BASES = ("A", "C", "G", "T")
MAX_SEED = 2**32 - 1


@dataclass(frozen=True)
class CandidateMutation:
    source: str
    pos: int
    new_base: str
    simulated_modern: str = "."
    simulated_ancient: str = "."


@dataclass(frozen=True)
class AppliedMutation:
    source: str
    haplotype: str
    pos: int
    original_base: str
    new_base: str
    simulated_modern: str
    simulated_ancient: str
    seed: int


def next_seed(seed: int, offset: int) -> int:
    return ((seed + offset - 1) % MAX_SEED) + 1


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


def mutation_model(name: str):
    if name == "jc69":
        return msprime.JC69()
    raise ValueError(f"Unsupported mutation model: {name}")


def branch_tree_sequence(sequence_length: int, generations: float) -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length=sequence_length)
    parent = tables.nodes.add_row(time=generations)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(left=0, right=sequence_length, parent=parent, child=child)
    tables.sort()
    return tables.tree_sequence()


def terminal_candidates(
    sequence_length: int,
    rate: float,
    generations: float,
    seed: int,
    model: str,
) -> list[CandidateMutation]:
    ts = branch_tree_sequence(sequence_length, generations)
    mts = msprime.sim_mutations(
        ts,
        rate=rate,
        model=mutation_model(model),
        random_seed=seed,
        discrete_genome=True,
    )

    out: list[CandidateMutation] = []
    sample_id = mts.samples()[0]
    for variant in mts.variants(samples=[sample_id]):
        pos = int(variant.site.position)
        allele = variant.alleles[int(variant.genotypes[0])]
        if allele in BASES:
            out.append(CandidateMutation("terminal_branch", pos, allele))
    return out


def split_candidates(
    sequence_length: int,
    rate: float,
    sample_generations: float,
    split_generations: float,
    modern_ne: float,
    ancient_ne: float,
    ancestral_ne: float,
    recombination_rate: float,
    seed: int,
    model: str,
) -> list[CandidateMutation]:
    if split_generations <= sample_generations:
        raise ValueError("--split-generations must be older than the ancient sample age")

    demography = msprime.Demography()
    demography.add_population(name="modern", initial_size=modern_ne)
    demography.add_population(name="ancient", initial_size=ancient_ne)
    demography.add_population(name="ancestral", initial_size=ancestral_ne)
    demography.add_population_split(
        time=split_generations,
        derived=["modern", "ancient"],
        ancestral="ancestral",
    )

    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(1, population="modern", time=0, ploidy=1),
            msprime.SampleSet(1, population="ancient", time=sample_generations, ploidy=1),
        ],
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )
    mts = msprime.sim_mutations(
        ts,
        rate=rate,
        model=mutation_model(model),
        random_seed=next_seed(seed, 1009),
        discrete_genome=True,
    )

    sample_nodes = list(mts.samples())
    if len(sample_nodes) != 2:
        raise ValueError(f"Expected 2 haploid samples from split model, found {len(sample_nodes)}")
    modern_node, ancient_node = sample_nodes

    out: list[CandidateMutation] = []
    for variant in mts.variants(samples=[modern_node, ancient_node]):
        pos = int(variant.site.position)
        modern_allele = variant.alleles[int(variant.genotypes[0])]
        ancient_allele = variant.alleles[int(variant.genotypes[1])]
        if modern_allele in BASES and ancient_allele in BASES and modern_allele != ancient_allele:
            out.append(
                CandidateMutation(
                    "split_demography",
                    pos,
                    ancient_allele,
                    simulated_modern=modern_allele,
                    simulated_ancient=ancient_allele,
                )
            )
    return out


def replacement_base(current: str, proposed: str, rng: random.Random) -> str:
    if proposed in BASES and proposed != current:
        return proposed
    choices = [base for base in BASES if base != current]
    return rng.choice(choices)


def diverge_sequence(
    header: str,
    sequence: str,
    candidates: list[CandidateMutation],
    seed: int,
) -> tuple[str, list[AppliedMutation], int]:
    rng = random.Random(seed + 104729)
    seq = list(sequence)
    applied: list[AppliedMutation] = []
    skipped_non_acgt = 0
    used_positions: set[int] = set()

    for candidate in sorted(candidates, key=lambda item: item.pos):
        if candidate.pos in used_positions:
            continue
        current = seq[candidate.pos].upper()
        if current not in BASES:
            skipped_non_acgt += 1
            continue
        new_base = replacement_base(current, candidate.new_base, rng)
        seq[candidate.pos] = new_base
        used_positions.add(candidate.pos)
        applied.append(
            AppliedMutation(
                candidate.source,
                header,
                candidate.pos,
                current,
                new_base,
                candidate.simulated_modern,
                candidate.simulated_ancient,
                seed,
            )
        )

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
    parser.add_argument("--mutation-report", type=Path)
    parser.add_argument("--header", default=None)
    parser.add_argument("--mode", default="terminal", choices=["terminal", "split"])
    parser.add_argument("--rate", required=True, type=float)
    parser.add_argument("--generations", required=True, type=float)
    parser.add_argument("--split-generations", type=float)
    parser.add_argument("--effective-pop-size", default=10000.0, type=float)
    parser.add_argument("--modern-ne", type=float)
    parser.add_argument("--ancient-ne", type=float)
    parser.add_argument("--ancestral-ne", type=float)
    parser.add_argument("--recombination-rate", default=1e-8, type=float)
    parser.add_argument("--model", default="jc69", choices=["jc69"])
    parser.add_argument("--seed", type=int)
    parser.add_argument("--pansn-delimiter", default="#")
    parser.add_argument("--require-pansn", action="store_true")
    args = parser.parse_args()

    if args.rate < 0:
        raise ValueError("--rate must be >= 0")
    if args.generations < 0:
        raise ValueError("--generations must be >= 0")
    modern_ne = args.modern_ne if args.modern_ne is not None else args.effective_pop_size
    ancient_ne = args.ancient_ne if args.ancient_ne is not None else args.effective_pop_size
    ancestral_ne = args.ancestral_ne if args.ancestral_ne is not None else args.effective_pop_size
    if args.effective_pop_size <= 0:
        raise ValueError("--effective-pop-size must be > 0")
    if modern_ne <= 0:
        raise ValueError("--modern-ne must be > 0")
    if ancient_ne <= 0:
        raise ValueError("--ancient-ne must be > 0")
    if ancestral_ne <= 0:
        raise ValueError("--ancestral-ne must be > 0")
    if args.recombination_rate < 0:
        raise ValueError("--recombination-rate must be >= 0")
    if args.mode == "split" and args.split_generations is None:
        raise ValueError("--mode split requires --split-generations")

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

    seed = args.seed if args.seed is not None else random.SystemRandom().randint(1, MAX_SEED)
    hap1_seed = seed
    hap2_seed = next_seed(seed, 1000003)

    if args.mode == "terminal":
        candidates1 = terminal_candidates(len(seq1), args.rate, args.generations, hap1_seed, args.model)
        candidates2 = terminal_candidates(len(seq2), args.rate, args.generations, hap2_seed, args.model)
    else:
        candidates1 = split_candidates(
            len(seq1),
            args.rate,
            args.generations,
            args.split_generations,
            modern_ne,
            ancient_ne,
            ancestral_ne,
            args.recombination_rate,
            hap1_seed,
            args.model,
        )
        candidates2 = split_candidates(
            len(seq2),
            args.rate,
            args.generations,
            args.split_generations,
            modern_ne,
            ancient_ne,
            ancestral_ne,
            args.recombination_rate,
            hap2_seed,
            args.model,
        )

    seq1_div, hap1_applied, hap1_skipped = diverge_sequence(header1, seq1, candidates1, hap1_seed)
    seq2_div, hap2_applied, hap2_skipped = diverge_sequence(header2, seq2, candidates2, hap2_seed)

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

    callable1 = sum(base in BASES for base in seq1)
    callable2 = sum(base in BASES for base in seq2)
    expected1 = args.rate * args.generations * callable1 if args.mode == "terminal" else "NA"
    expected2 = args.rate * args.generations * callable2 if args.mode == "terminal" else "NA"

    if args.report is not None:
        with args.report.open("w") as handle:
            handle.write(
                "haplotype\tmode\tsource\tlength\tcallable_acgt\texpected_mutations\t"
                "candidate_mutations\tapplied_mutations\tskipped_non_acgt\tseed\t"
                "sample_generations\tsplit_generations\tmodern_ne\tancient_ne\t"
                "ancestral_ne\trecombination_rate\n"
            )
            for header, seq, callable_bases, expected, candidates, applied, skipped, row_seed in (
                (header1, seq1, callable1, expected1, candidates1, hap1_applied, hap1_skipped, hap1_seed),
                (header2, seq2, callable2, expected2, candidates2, hap2_applied, hap2_skipped, hap2_seed),
            ):
                source = "terminal_branch" if args.mode == "terminal" else "split_demography"
                handle.write(
                    f"{header}\t{args.mode}\t{source}\t{len(seq)}\t{callable_bases}\t"
                    f"{expected}\t{len(candidates)}\t{len(applied)}\t{skipped}\t{row_seed}\t"
                    f"{args.generations}\t{args.split_generations if args.split_generations is not None else 'NA'}\t"
                    f"{modern_ne}\t{ancient_ne}\t{ancestral_ne}\t{args.recombination_rate}\n"
                )

    if args.mutation_report is not None:
        with args.mutation_report.open("w") as handle:
            handle.write(
                "source\thaplotype\tposition_0based\tposition_1based\toriginal_base\tnew_base\t"
                "simulated_modern\tsimulated_ancient\tseed\n"
            )
            for item in hap1_applied + hap2_applied:
                handle.write(
                    f"{item.source}\t{item.haplotype}\t{item.pos}\t{item.pos + 1}\t"
                    f"{item.original_base}\t{item.new_base}\t{item.simulated_modern}\t"
                    f"{item.simulated_ancient}\t{item.seed}\n"
                )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
