"""Sequence-based paralog family detection from the UniProt N. crassa proteome.

Novel extension using data/raw/uniprotkb_proteome_UP000001805_2026_07_02.fasta
(the UniProt reference proteome for N. crassa, UP000001805 -- 10,255 protein
sequences). This is NOT part of the original paper; it answers a different
question the paper doesn't address at all: which clock-network genes belong
to small (size-3 or size-4) paralogous gene families? That's the kind of
family definition needed to apply a duplication-aware framework (the
companion Scruse/Arnold/Robinson Partial-Duplication project) to this
network -- something not possible without real sequence data.

Method, in two stages (no BLAST/diamond/mmseqs/cd-hit available in this
environment):

1. **Cheap prefilter**: protein k-mer (default k=5) Jaccard similarity via
   an inverted-index, generating candidate pairs. This stage is deliberately
   loose (low jaccard threshold) -- its only job is to shrink ~10,255^2
   possible pairs down to a few hundred plausible candidates, cheaply.
2. **Real confirmation**: each candidate pair is then actually aligned with
   Biopython's PairwiseAligner (local alignment, BLOSUM62, affine gaps) --
   the same substitution-matrix-based approach BLAST itself is built on,
   just without BLAST's statistical e-value machinery. A pair is called
   "same family" only if percent identity and alignment coverage both clear
   a threshold on this real alignment, not on the k-mer heuristic.

This two-stage design exists because stage 1 alone is a poor family
definition in practice: at a strict k-mer threshold it mostly recovers
duplicate/redundant UniProt records of the *same* gene (not paralogs), and at
a loose threshold it's dominated by spurious coincidental k-mer overlaps.
Real alignment in stage 2 is what actually decides family membership.

Known limitation: a k-mer prefilter can still miss true but highly divergent
(remote) paralogs that happen to share very few exact k-mers -- this is not
a substitute for a profile/HMM-based search (e.g. OrthoFinder, InterPro).
Single-linkage clustering of the surviving edges can also over-merge distinct
families through a promiscuous intermediate. Treat results as a first-pass,
not a definitive orthology/paralogy call.

Identifier resolution caveat: this FASTA uses a mix of NCU-locus-tag and
common gene-symbol GN= values, and cross-checking against the paper's own
regulator names surfaced real discrepancies -- see REGULATOR_ALIASES and
UNRESOLVED_REGULATORS below, both determined by direct inspection of the
file, not assumed.
"""
from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations

HEADER_RE = re.compile(r"^>(\w+)\|([^|]+)\|(\S+)\s+(.*)$")
GENE_RE = re.compile(r"\bGN=(\S+)")
DESC_SPLIT_RE = re.compile(r"\s(?:OS|OX|GN|PE|SV)=")
NCU_PATTERN = re.compile(r"^NCU\d{5}$")

# Resolved by direct inspection of the FASTA: the paper's regulator common
# name doesn't match UniProt's GN= value, but the UniProt entry's function
# description confirms it's the same gene.
REGULATOR_ALIASES = {
    "NCU04799": "pabp-1",  # paper calls this "pab-1"; UniProt entry PABP_NEUCR, GN=pabp-1 -- same protein (polyadenylate-binding protein)
    "NCU00919": "drh-16",  # paper calls this "rok-1"; UniProt entry ROK1_NEUCR ("ATP-dependent RNA helicase rok1"), GN=drh-16 -- same protein, different symbol on record
}

# Regulators that could NOT be confidently resolved to a sequence in this
# FASTA -- excluded from family analysis rather than guessing.
UNRESOLVED_REGULATORS = {
    "NCU03363": (
        "No matching UniProt entry found by NCU locus tag or any name variant "
        "tried (searched the raw file directly). Likely a locus retired/renumbered "
        "in a genome reannotation between this paper (2022) and the current "
        "UniProt reference proteome."
    ),
    "NCU07392": (
        "The paper labels this regulator \"pro-1\", but UniProt's GN=pro-1 "
        "(Q12641) is annotated as Pyrroline-5-carboxylate reductase -- a "
        "well-established, distinct classical Neurospora proline-biosynthesis "
        "locus, not a DNA-binding transcriptional regulator. This looks like a "
        "gene-symbol collision or a naming error rather than the same gene; "
        "excluded here rather than risk building families from the wrong "
        "sequence. Needs independent verification against the paper's original "
        "gene models to resolve."
    ),
}

# WCC is the WC-1/WC-2 protein complex, not a single gene product -- no
# single sequence to assign it, so it's excluded from family analysis (its
# two subunits, wc-1 and wc-2, are core clock genes handled separately by
# data_loader.py and are not part of the target-gene universe either).
EXCLUDED_REGULATORS = {"Wcc": "WCC is the WC-1/WC-2 protein complex, not a single gene product."}


@dataclass
class ProteinRecord:
    uniprot_id: str
    entry_name: str
    gene_name: str | None
    description: str
    sequence: str


@dataclass
class FamilyBuildReport:
    n_records_parsed: int
    n_records_deduped: int
    n_stage1_candidates: int
    n_stage2_confirmed_edges: int
    n_duplicate_record_edges: int
    n_paralog_edges: int
    n_families: int
    family_size_counts: dict  # size -> number of families of that size
    kmer_size: int
    stage1_jaccard_threshold: float
    stage1_min_shared_kmers: int
    stage1_max_posting_list: int
    stage2_min_identity: float
    stage2_min_coverage: float


def parse_fasta(path: str) -> list[ProteinRecord]:
    records: list[ProteinRecord] = []
    header = None
    seq_lines: list[str] = []

    def flush():
        if header is None:
            return
        db, uid, entry, rest = header
        m = DESC_SPLIT_RE.search(rest)
        desc = rest[: m.start()] if m else rest
        gene_m = GENE_RE.search(rest)
        gene = gene_m.group(1) if gene_m else None
        records.append(ProteinRecord(uid, entry, gene, desc.strip(), "".join(seq_lines)))

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                m = HEADER_RE.match(line)
                if not m:
                    header = None
                    seq_lines = []
                    continue
                header = (m.group(1), m.group(2), m.group(3), m.group(4))
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        flush()
    return records


def dedupe_by_gene_name(records: list[ProteinRecord]) -> dict[str, ProteinRecord]:
    """Keep one record per gene name (the longest sequence, if duplicated)."""
    best: dict[str, ProteinRecord] = {}
    for r in records:
        if r.gene_name is None:
            continue
        prev = best.get(r.gene_name)
        if prev is None or len(r.sequence) > len(prev.sequence):
            best[r.gene_name] = r
    return best


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def _kmer_set(seq: str, k: int) -> frozenset:
    if len(seq) < k:
        return frozenset()
    return frozenset(seq[i : i + k] for i in range(len(seq) - k + 1))


def stage1_kmer_candidates(
    records: dict[str, ProteinRecord],
    k: int = 5,
    max_posting_list: int = 40,
    min_shared_kmers: int = 3,
    jaccard_threshold: float = 0.05,
) -> list[tuple[str, str, float]]:
    """Cheap prefilter: candidate pairs by protein k-mer Jaccard similarity.

    Deliberately loose -- only meant to shrink the search space before real
    alignment in stage 2, not to make the final family call.
    """
    gene_names = list(records.keys())
    kmer_sets = {g: _kmer_set(records[g].sequence, k) for g in gene_names}

    inverted: dict[str, list[str]] = defaultdict(list)
    for gene, kmers in kmer_sets.items():
        for km in kmers:
            inverted[km].append(gene)

    candidate_pairs: set[tuple[str, str]] = set()
    for gene_list in inverted.values():
        if len(gene_list) < 2 or len(gene_list) > max_posting_list:
            continue
        uniq = sorted(set(gene_list))
        for a, b in combinations(uniq, 2):
            candidate_pairs.add((a, b))

    results = []
    for a, b in candidate_pairs:
        A, B = kmer_sets[a], kmer_sets[b]
        inter = len(A & B)
        if inter < min_shared_kmers:
            continue
        union_size = len(A) + len(B) - inter
        if union_size == 0:
            continue
        jac = inter / union_size
        if jac >= jaccard_threshold:
            results.append((a, b, jac))
    return results


def stage2_alignment_confirm(
    candidates: list[tuple[str, str, float]],
    records: dict[str, ProteinRecord],
    min_identity: float = 0.25,
    min_coverage: float = 0.5,
) -> list[tuple[str, str, float, float]]:
    """Confirm stage-1 candidates with real local alignment (BLOSUM62).

    Returns (gene_a, gene_b, percent_identity, coverage) for pairs clearing
    both thresholds, where coverage = aligned length / shorter sequence's
    length (guards against a short shared motif/domain inflating identity%
    on an otherwise unrelated pair of proteins).
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    confirmed = []
    for a, b, _jac in candidates:
        seq_a, seq_b = records[a].sequence, records[b].sequence
        if not seq_a or not seq_b:
            continue
        aln = aligner.align(seq_a, seq_b)[0]
        counts = aln.counts()
        if counts.aligned == 0:
            continue
        identity = counts.identities / counts.aligned
        coverage = counts.aligned / min(len(seq_a), len(seq_b))
        if identity >= min_identity and coverage >= min_coverage:
            confirmed.append((a, b, identity, coverage))
    return confirmed


def cluster_families(gene_names: list[str], edges: list[tuple]) -> dict[str, int]:
    """Union-find single-linkage clustering. Returns gene_name -> family_id."""
    uf = UnionFind(gene_names)
    for edge in edges:
        a, b = edge[0], edge[1]
        uf.union(a, b)

    roots: dict[str, int] = {}
    family_of: dict[str, int] = {}
    next_id = 0
    for g in gene_names:
        r = uf.find(g)
        if r not in roots:
            roots[r] = next_id
            next_id += 1
        family_of[g] = roots[r]
    return family_of


# A confirmed edge at (near-)100% identity and coverage is almost always the
# *same physical gene* recorded twice under two different UniProt
# accessions (one NCU-locus-tag entry, one common-name entry) -- not a real
# paralog pair. Observed directly in this proteome: e.g. matA-1 <-> mtA-1
# and NCU02680 <-> npl4 both align at ~100% identity/coverage and are the
# same annotated function under two names. These must be collapsed to a
# single canonical gene *before* family clustering, or "family size" would
# be inflated by annotation duplication rather than real duplication.
DUPLICATE_RECORD_IDENTITY = 0.97
DUPLICATE_RECORD_COVERAGE = 0.85


def split_duplicate_records(
    confirmed_edges: list[tuple[str, str, float, float]],
) -> tuple[list[tuple], list[tuple]]:
    """Split confirmed edges into (duplicate-record edges, genuine paralog edges)."""
    dup_edges, paralog_edges = [], []
    for edge in confirmed_edges:
        _a, _b, identity, coverage = edge
        if identity >= DUPLICATE_RECORD_IDENTITY and coverage >= DUPLICATE_RECORD_COVERAGE:
            dup_edges.append(edge)
        else:
            paralog_edges.append(edge)
    return dup_edges, paralog_edges


def canonicalize_duplicate_records(
    gene_names: list[str], records: dict[str, ProteinRecord], dup_edges: list[tuple]
) -> tuple[dict[str, str], dict[str, ProteinRecord]]:
    """Collapse duplicate-record groups into one canonical gene name each.

    Prefers an NCU-locus-tag name as canonical (so downstream lookups by
    NCU id keep working); otherwise keeps the first name encountered.
    Returns (original_gene_name -> canonical_name, canonical_name -> record).
    """
    uf = UnionFind(gene_names)
    for a, b, _i, _c in dup_edges:
        uf.union(a, b)

    groups: dict[str, list[str]] = defaultdict(list)
    for g in gene_names:
        groups[uf.find(g)].append(g)

    canonical_of: dict[str, str] = {}
    canonical_records: dict[str, ProteinRecord] = {}
    for members in groups.values():
        ncu_members = sorted(m for m in members if NCU_PATTERN.match(m))
        canonical = ncu_members[0] if ncu_members else sorted(members)[0]
        best_record = max((records[m] for m in members), key=lambda r: len(r.sequence))
        for m in members:
            canonical_of[m] = canonical
        canonical_records[canonical] = best_record
    return canonical_of, canonical_records


def build_families(
    records: dict[str, ProteinRecord],
    k: int = 5,
    max_posting_list: int = 40,
    stage1_min_shared_kmers: int = 3,
    stage1_jaccard_threshold: float = 0.05,
    stage2_min_identity: float = 0.25,
    stage2_min_coverage: float = 0.5,
) -> tuple[dict[str, int], dict[str, str], list[tuple[str, str, float, float]], FamilyBuildReport]:
    """Full pipeline: k-mer prefilter -> real alignment confirmation ->
    collapse same-gene duplicate UniProt records -> single-linkage cluster
    the *remaining* genuine-paralog edges into families.

    Returns (canonical_gene_name -> family_id, original_name -> canonical_name,
    paralog_edges [on canonical names], report).
    """
    candidates = stage1_kmer_candidates(
        records, k=k, max_posting_list=max_posting_list,
        min_shared_kmers=stage1_min_shared_kmers, jaccard_threshold=stage1_jaccard_threshold,
    )
    confirmed = stage2_alignment_confirm(
        candidates, records, min_identity=stage2_min_identity, min_coverage=stage2_min_coverage,
    )
    dup_edges, paralog_edges = split_duplicate_records(confirmed)
    canonical_of, canonical_records = canonicalize_duplicate_records(
        list(records.keys()), records, dup_edges
    )

    remapped_paralog_edges = []
    for a, b, identity, coverage in paralog_edges:
        ca, cb = canonical_of.get(a, a), canonical_of.get(b, b)
        if ca != cb:
            remapped_paralog_edges.append((ca, cb, identity, coverage))

    canonical_names = list(canonical_records.keys())
    family_of = cluster_families(canonical_names, remapped_paralog_edges)

    family_sizes: dict[int, int] = defaultdict(int)
    for fam_id in family_of.values():
        family_sizes[fam_id] += 1
    size_counts: dict[int, int] = defaultdict(int)
    for sz in family_sizes.values():
        size_counts[sz] += 1

    report = FamilyBuildReport(
        n_records_parsed=len(records),
        n_records_deduped=len(canonical_names),
        n_stage1_candidates=len(candidates),
        n_stage2_confirmed_edges=len(confirmed),
        n_duplicate_record_edges=len(dup_edges),
        n_paralog_edges=len(remapped_paralog_edges),
        n_families=len(family_sizes),
        family_size_counts=dict(sorted(size_counts.items())),
        kmer_size=k,
        stage1_jaccard_threshold=stage1_jaccard_threshold,
        stage1_min_shared_kmers=stage1_min_shared_kmers,
        stage1_max_posting_list=max_posting_list,
        stage2_min_identity=stage2_min_identity,
        stage2_min_coverage=stage2_min_coverage,
    )
    return family_of, canonical_of, remapped_paralog_edges, report


def resolve_regulator_gene_names(regulator_info: dict) -> dict[str, str | None]:
    """Map each regulator label (e.g. "PostReg9 NCU00919") to the gene-name
    key it should be looked up under in the parsed FASTA records, applying
    the alias table above. Returns None for regulators that are excluded or
    unresolved (see UNRESOLVED_REGULATORS / EXCLUDED_REGULATORS)."""
    resolved = {}
    for label, info in regulator_info.items():
        ncu = info.get("ncu")
        if label in EXCLUDED_REGULATORS:
            resolved[label] = None
        elif ncu in UNRESOLVED_REGULATORS:
            resolved[label] = None
        elif ncu in REGULATOR_ALIASES:
            resolved[label] = REGULATOR_ALIASES[ncu]
        elif ncu is not None:
            resolved[label] = ncu
        else:
            resolved[label] = None
    return resolved


DEFAULT_FASTA_PATH = "data/raw/uniprotkb_proteome_UP000001805_2026_07_02.fasta"
DEFAULT_OUTPUT_PATH = "data/processed/gene_families.csv"


def summarize_k34_families(
    family_of: dict[str, int],
    canonical_records: dict[str, ProteinRecord],
    clock_gene_ids: set[str],
) -> str:
    groups: dict[int, list[str]] = defaultdict(list)
    for g, fid in family_of.items():
        groups[fid].append(g)

    lines = []
    for k in (3, 4):
        fams = [members for members in groups.values() if len(members) == k]
        lines.append(f"\n--- Families of size k={k} ({len(fams)} found) ---")
        for members in fams:
            touches_clock = [m for m in members if m in clock_gene_ids]
            marker = f"  <-- includes clock-network gene(s): {touches_clock}" if touches_clock else ""
            lines.append(f"  {members}{marker}")
            for m in members:
                lines.append(f"      {m}: {canonical_records[m].description}")
    return "\n".join(lines)


if __name__ == "__main__":
    import os
    import time

    import pandas as pd

    from data_loader import REGULATOR_INFO
    from network_builder import build_network

    t0 = time.time()
    records = parse_fasta(DEFAULT_FASTA_PATH)
    print(f"Parsed {len(records)} records in {time.time() - t0:.1f}s")

    deduped = dedupe_by_gene_name(records)
    print(f"Deduped to {len(deduped)} distinct gene names")

    t0 = time.time()
    family_of, canonical_of, paralog_edges, report = build_families(deduped)
    print(f"Built families in {time.time() - t0:.1f}s")
    print(report)

    canonical_records = {g: deduped[g] for g in family_of}

    print(f"\nGenuine paralog edges used for family clustering ({len(paralog_edges)}):")
    for a, b, ident, cov in sorted(paralog_edges, key=lambda e: -e[2]):
        print(f"  {a} <-> {b}: identity={ident:.2f} coverage={cov:.2f} "
              f"[{canonical_records[a].description} | {canonical_records[b].description}]")

    # Cross-reference against the clock network.
    net = build_network()
    regulator_gene_names = resolve_regulator_gene_names(REGULATOR_INFO)
    clock_gene_ids = set(net.T.columns)
    for label, gene_name in regulator_gene_names.items():
        if gene_name:
            canonical = canonical_of.get(gene_name, gene_name)
            clock_gene_ids.add(canonical)

    print(summarize_k34_families(family_of, canonical_records, clock_gene_ids))

    print("\n--- Regulator resolution status ---")
    for label, gene_name in regulator_gene_names.items():
        ncu = REGULATOR_INFO[label]["ncu"]
        paper_name = REGULATOR_INFO[label]["gene"]
        if gene_name is None:
            reason = (
                UNRESOLVED_REGULATORS.get(ncu)
                or EXCLUDED_REGULATORS.get(label)
                or "unresolved"
            )
            print(f"  {label} ({paper_name}): NOT RESOLVED -- {reason}")
        else:
            canonical = canonical_of.get(gene_name, gene_name)
            fam_id = family_of.get(canonical)
            fam_size = sum(1 for v in family_of.values() if v == fam_id) if fam_id is not None else 1
            print(f"  {label} ({paper_name}) -> {gene_name} -> family size {fam_size}")

    # Cache full gene -> family assignment for the app.
    os.makedirs(os.path.dirname(DEFAULT_OUTPUT_PATH), exist_ok=True)
    family_sizes = defaultdict(int)
    for fid in family_of.values():
        family_sizes[fid] += 1
    out_df = pd.DataFrame(
        {
            "gene_name": list(family_of.keys()),
            "family_id": list(family_of.values()),
            "description": [canonical_records[g].description for g in family_of],
            "is_clock_network_gene": [g in clock_gene_ids for g in family_of],
        }
    )
    out_df["family_size"] = out_df["family_id"].map(family_sizes)
    out_df = out_df.sort_values(["family_size", "family_id"], ascending=[False, True])
    out_df.to_csv(DEFAULT_OUTPUT_PATH, index=False)
    print(f"\nSaved gene->family assignments to {DEFAULT_OUTPUT_PATH}")
