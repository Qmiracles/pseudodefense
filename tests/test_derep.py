import gzip
from pathlib import Path

import pandas as pd
import pytest

from pseudodefense import derep

# Distances (1 - ANI): A-B 0.0005 (99.95 % ANI), C 0.008 from A/B (99.2 %), D 0.02 from all
PAIRS = {
    ("A", "B"): 0.0005,
    ("A", "C"): 0.008,
    ("B", "C"): 0.008,
    ("A", "D"): 0.02,
    ("B", "D"): 0.02,
    ("C", "D"): 0.02,
}


def _write_mash(path: Path) -> Path:
    names = ["A", "B", "C", "D"]
    lines = []
    for r in names:
        for q in names:
            d = 0.0 if r == q else PAIRS.get((r, q), PAIRS.get((q, r)))
            lines.append(f"results/genomes/{r}.fna.gz\tresults/genomes/{q}.fna.gz\t{d}\t0\t9/10\n")
    with gzip.open(path, "wt") as fh:
        fh.writelines(lines)
    return path


def test_accession_from_path():
    assert derep.accession_from_path("x/GCF_000006765.1.fna.gz") == "GCF_000006765.1"


def test_read_mash_dist_square_and_symmetric(tmp_path):
    d = derep.read_mash_dist(_write_mash(tmp_path / "d.tsv.gz"))
    assert d.index.tolist() == ["A", "B", "C", "D"]
    assert d.loc["A", "B"] == pytest.approx(0.0005)
    assert (d.to_numpy() == d.to_numpy().T).all()


def test_read_mash_dist_incomplete_fails(tmp_path):
    path = tmp_path / "d.tsv"
    path.write_text("A.fna.gz\tB.fna.gz\t0.1\t0\t1/10\nA.fna.gz\tA.fna.gz\t0\t0\t10/10\n")
    with pytest.raises(ValueError, match="incomplete"):
        derep.read_mash_dist(path)


def test_cluster_genomes_by_threshold(tmp_path):
    d = derep.read_mash_dist(_write_mash(tmp_path / "d.tsv.gz"))
    clusters = derep.cluster_genomes(d, [99.9, 99.0]).set_index("accession")
    # 99.9 %: only A and B merge
    assert clusters.loc["A", "ani99.9"] == clusters.loc["B", "ani99.9"]
    assert clusters["ani99.9"].nunique() == 3
    # 99.0 %: A, B and C merge; D stays alone
    assert clusters.loc[["A", "B", "C"], "ani99"].nunique() == 1
    assert clusters["ani99"].nunique() == 2


def test_representatives_prefer_quality():
    clusters = pd.DataFrame({"accession": ["A", "B", "C"], "ani99": [1, 1, 2]})
    manifest = pd.DataFrame(
        {
            "accession": ["A", "B", "C"],
            "assembly_level": ["Chromosome", "Complete Genome", "Chromosome"],
        }
    )
    qc_report = pd.DataFrame(
        {
            "accession": ["A", "B", "C"],
            "n_sequences": [1, 3, 1],
            "checkm_completeness": [99.9, 98.0, None],
            "checkm_contamination": [0.1, 0.5, None],
        }
    )
    ranks = derep.quality_order(manifest, qc_report)
    reps = derep.choose_representatives(clusters, ranks, [99.0]).set_index("accession")
    # B is Complete Genome, so it wins over A despite lower CheckM
    assert reps["rep_ani99"].to_dict() == {"A": False, "B": True, "C": True}


def test_summary_counts():
    clusters = pd.DataFrame({"accession": list("ABCD"), "ani99": [1, 1, 1, 2]})
    s = derep.summary(clusters, [99.0]).iloc[0]
    assert (s.n_clusters, s.largest_cluster, s.singletons) == (2, 3, 1)


def test_select_validation_pairs_near_thresholds(tmp_path):
    d = derep.read_mash_dist(_write_mash(tmp_path / "d.tsv.gz"))
    pairs = derep.select_validation_pairs(d, [99.95, 99.2], per_threshold=5, window=0.01, seed=1)
    got = {(r.accession_a, r.accession_b, r.near_threshold) for r in pairs.itertuples()}
    # A-B sits at 99.95 %; A-C and B-C at 99.2 %; nothing else is inside the bands
    assert got == {("A", "B", 99.95), ("A", "C", 99.2), ("B", "C", 99.2)}


def test_select_validation_pairs_respects_sample_size(tmp_path):
    d = derep.read_mash_dist(_write_mash(tmp_path / "d.tsv.gz"))
    pairs = derep.select_validation_pairs(d, [99.2], per_threshold=1, window=0.01, seed=1)
    assert len(pairs) == 1


def test_read_fastani_and_compare(tmp_path):
    raw = tmp_path / "fastani.tsv"
    # FastANI may report the pair in either order; C-A falls on the other side of 99.2
    raw.write_text(
        "g/B.fna.gz\tg/A.fna.gz\t99.96\t2000\t2100\ng/C.fna.gz\tg/A.fna.gz\t99.1\t9\t9\n"
    )
    pairs = pd.DataFrame(
        {
            "accession_a": ["A", "A", "B"],
            "accession_b": ["B", "C", "C"],
            "mash_ani": [99.95, 99.2, 99.2],
            "near_threshold": [99.95, 99.2, 99.2],
        }
    )
    table, summ = derep.compare_ani(pairs, derep.read_fastani(raw))
    assert table.loc[0, "difference"] == pytest.approx(0.01)
    assert pd.isna(table.loc[2, "fastani"])  # B-C missing from FastANI output
    s = summ.set_index("threshold_ani")
    assert s.loc[99.2, "n_with_fastani"] == 1
    assert s.loc[99.2, "same_side_fraction"] == 0.0
    assert s.loc[99.95, "same_side_fraction"] == 1.0


def test_read_fastani_empty(tmp_path):
    raw = tmp_path / "fastani.tsv"
    raw.write_text("")
    assert derep.read_fastani(raw).empty
