"""Tests for sequence_families.py, using small synthetic inputs (not the full
proteome -- the real pipeline takes ~3 minutes and is exercised manually via
`python src/sequence_families.py`, not in the test suite)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sequence_families import (  # noqa: E402
    ProteinRecord,
    UnionFind,
    build_families,
    canonicalize_duplicate_records,
    cluster_families,
    dedupe_by_gene_name,
    parse_fasta,
    resolve_regulator_gene_names,
    split_duplicate_records,
    stage1_kmer_candidates,
    stage2_alignment_confirm,
)

FASTA_FIXTURE = """\
>sp|P00001|GENEA_NEUCR Some enzyme OS=Neurospora crassa OX=367110 GN=gene-a PE=1 SV=1
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPR
GRRQPIPKAR
>sp|P00002|GENEB_NEUCR Some enzyme, duplicate record OS=Neurospora crassa OX=367110 GN=NCU00002 PE=1 SV=1
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKAR
>sp|P00003|GENEC_NEUCR Unrelated protein OS=Neurospora crassa OX=367110 GN=NCU00003 PE=1 SV=1
GGGGVVVVLLLLIIIIAAAAWWWWFFFFYYYYMMMMCCCCPPPPSSSSTTTTNNNNQQQ
>sp|P00004|GENED_NEUCR No gene name here OS=Neurospora crassa OX=367110 PE=1 SV=1
ACDEFGHIKLMNPQRSTVWY
"""


@pytest.fixture
def fasta_path(tmp_path):
    p = tmp_path / "mini.fasta"
    p.write_text(FASTA_FIXTURE, encoding="utf-8")
    return str(p)


def test_parse_fasta_extracts_gene_names_and_sequences(fasta_path):
    records = parse_fasta(fasta_path)
    assert len(records) == 4
    by_uid = {r.uniprot_id: r for r in records}

    assert by_uid["P00001"].gene_name == "gene-a"
    assert by_uid["P00001"].sequence == "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKAR"
    assert by_uid["P00001"].description == "Some enzyme"

    assert by_uid["P00002"].gene_name == "NCU00002"
    assert by_uid["P00004"].gene_name is None


def test_dedupe_by_gene_name_drops_records_without_gene_name(fasta_path):
    records = parse_fasta(fasta_path)
    deduped = dedupe_by_gene_name(records)
    assert set(deduped.keys()) == {"gene-a", "NCU00002", "NCU00003"}


def test_dedupe_keeps_longer_sequence_on_collision():
    records = [
        ProteinRecord("A1", "A1_X", "dup-gene", "desc", "SHORT"),
        ProteinRecord("A2", "A2_X", "dup-gene", "desc", "MUCHLONGERSEQUENCE"),
    ]
    deduped = dedupe_by_gene_name(records)
    assert deduped["dup-gene"].sequence == "MUCHLONGERSEQUENCE"


def test_stage1_finds_near_identical_pair_but_not_unrelated_one(fasta_path):
    records = dedupe_by_gene_name(parse_fasta(fasta_path))
    candidates = stage1_kmer_candidates(records, k=5, min_shared_kmers=3, jaccard_threshold=0.05)
    pairs = {frozenset((a, b)) for a, b, _j in candidates}
    assert frozenset(("gene-a", "NCU00002")) in pairs
    assert frozenset(("gene-a", "NCU00003")) not in pairs


def test_stage2_confirms_identical_sequences_at_full_identity(fasta_path):
    records = dedupe_by_gene_name(parse_fasta(fasta_path))
    candidates = [("gene-a", "NCU00002", 1.0)]
    confirmed = stage2_alignment_confirm(candidates, records, min_identity=0.9, min_coverage=0.9)
    assert len(confirmed) == 1
    a, b, identity, coverage = confirmed[0]
    assert identity == pytest.approx(1.0, abs=0.01)
    assert coverage == pytest.approx(1.0, abs=0.05)


def test_stage2_rejects_unrelated_sequences(fasta_path):
    records = dedupe_by_gene_name(parse_fasta(fasta_path))
    candidates = [("gene-a", "NCU00003", 0.5)]  # forced candidate despite no real similarity
    confirmed = stage2_alignment_confirm(candidates, records, min_identity=0.25, min_coverage=0.5)
    assert confirmed == []


def test_union_find_basic():
    uf = UnionFind(["a", "b", "c", "d"])
    uf.union("a", "b")
    uf.union("b", "c")
    assert uf.find("a") == uf.find("c")
    assert uf.find("a") != uf.find("d")


def test_cluster_families_groups_transitively():
    genes = ["a", "b", "c", "d"]
    edges = [("a", "b", 0.9, 0.9), ("b", "c", 0.8, 0.8)]
    family_of = cluster_families(genes, edges)
    assert family_of["a"] == family_of["b"] == family_of["c"]
    assert family_of["d"] != family_of["a"]


def test_split_duplicate_records_separates_by_identity_and_coverage():
    edges = [
        ("a", "b", 0.99, 0.98),  # duplicate record
        ("c", "d", 0.60, 0.90),  # genuine paralog
        ("e", "f", 0.99, 0.50),  # high identity but low coverage -> not a duplicate record
    ]
    dup_edges, paralog_edges = split_duplicate_records(edges)
    assert ("a", "b", 0.99, 0.98) in dup_edges
    assert ("c", "d", 0.60, 0.90) in paralog_edges
    assert ("e", "f", 0.99, 0.50) in paralog_edges


def test_canonicalize_prefers_ncu_pattern_name():
    records = {
        "friendly-name": ProteinRecord("U1", "E1", "friendly-name", "d", "AAAA"),
        "NCU00099": ProteinRecord("U2", "E2", "NCU00099", "d", "AAAAA"),
    }
    dup_edges = [("friendly-name", "NCU00099", 1.0, 1.0)]
    canonical_of, canonical_records = canonicalize_duplicate_records(
        list(records.keys()), records, dup_edges
    )
    assert canonical_of["friendly-name"] == "NCU00099"
    assert canonical_of["NCU00099"] == "NCU00099"
    assert set(canonical_records.keys()) == {"NCU00099"}


def test_build_families_end_to_end_on_synthetic_set():
    # 3 near-identical "paralogs" (family of 3), 1 unrelated singleton, and a
    # same-gene duplicate-record pair that must NOT count as a real family.
    base = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKAR"
    # ~88% identical to base (scattered substitutions every ~9 residues) --
    # diverged enough to be a genuine paralog (well under the 97% duplicate-
    # record cutoff) while still leaving long enough conserved runs to share
    # real 5-mers with base, so the stage-1 k-mer prefilter actually finds it.
    mutated = "MSDNPKPQRKTERNTNRRPQFVKFPGGGQGVGGVYLLPHRGPRLGVRITRKTSERSKPRGRRQPILKAR"
    unrelated = "GGGGVVVVLLLLIIIIAAAAWWWWFFFFYYYYMMMMCCCCPPPPSSSSTTTTNNNNQQQQEEEEDDDD"

    records = {
        "para-1": ProteinRecord("U1", "E1", "para-1", "d1", base),
        "para-2": ProteinRecord("U2", "E2", "para-2", "d2", mutated),
        "NCU00001": ProteinRecord("U3", "E3", "NCU00001", "d3", base),  # exact duplicate of para-1 -> same gene
        "NCU00002": ProteinRecord("U4", "E4", "NCU00002", "d4", unrelated),
    }
    family_of, canonical_of, paralog_edges, report = build_families(records, k=5)

    # para-1 and NCU00001 are the *same* gene (100% identical) -> collapsed to one canonical id
    assert canonical_of["para-1"] == canonical_of["NCU00001"]
    canonical_para1 = canonical_of["para-1"]

    # that canonical gene and para-2 are genuine (imperfect) paralogs -> same family
    assert family_of[canonical_para1] == family_of["para-2"]

    # the unrelated gene must not join that family
    assert family_of["NCU00002"] != family_of[canonical_para1]

    assert report.n_duplicate_record_edges >= 1


def test_resolve_regulator_gene_names_handles_alias_unresolved_and_excluded():
    info = {
        "Wcc": {"gene": "WCC", "ncu": None, "class": "master"},
        "Reg1 NCU07392": {"gene": "pro-1", "ncu": "NCU07392", "class": "transcriptional"},
        "PostReg8 NCU04799": {"gene": "pab-1", "ncu": "NCU04799", "class": "post-transcriptional"},
        "PostReg5 NCU08295": {"gene": "lhp-1", "ncu": "NCU08295", "class": "post-transcriptional"},
    }
    resolved = resolve_regulator_gene_names(info)
    assert resolved["Wcc"] is None  # excluded, protein complex
    assert resolved["Reg1 NCU07392"] is None  # unresolved, flagged naming collision
    assert resolved["PostReg8 NCU04799"] == "pabp-1"  # alias
    assert resolved["PostReg5 NCU08295"] == "NCU08295"  # direct NCU match
