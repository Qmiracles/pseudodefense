import numpy as np
import pandas as pd
import pytest

from pseudodefense import stats


def _matrix(columns: dict[str, list[int]]) -> pd.DataFrame:
    n = len(next(iter(columns.values())))
    return pd.DataFrame(columns, index=[f"G{i}" for i in range(n)])


def test_wilson_interval_zero_successes():
    # Upper bound for 0/20: z^2/n / (1 + z^2/n) = 0.1611
    low, high = stats.wilson_interval(0, 20)
    assert low == pytest.approx(0)
    assert high == pytest.approx(0.1611, abs=1e-4)


def test_wilson_interval_half():
    low, high = stats.wilson_interval(10, 20)
    assert low == pytest.approx(1 - high)
    assert low == pytest.approx(0.2993, abs=1e-4)


def test_prevalence_sorted_and_bounded():
    m = _matrix({"A": [1, 1, 1, 0], "B": [1, 0, 0, 0], "C": [0, 0, 0, 0]})
    prev = stats.prevalence(m)
    assert prev["system"].tolist() == ["A", "B", "C"]
    assert prev["prevalence"].tolist() == [0.75, 0.25, 0.0]
    assert (prev["ci95_low"] <= prev["prevalence"]).all()
    assert (prev["ci95_high"] >= prev["prevalence"]).all()


def test_testable_systems_excludes_core_and_absent():
    m = _matrix({"core": [1] * 10, "absent": [0] * 10, "mid": [1] * 5 + [0] * 5})
    assert stats.testable_systems(m, 0.05, 0.95) == ["mid"]


def test_testable_systems_min_genomes():
    # 1/20 = 5 % passes the prevalence filter but not min_genomes=3
    m = _matrix({"rare": [1] + [0] * 19, "mid": [1] * 3 + [0] * 17, "common": [1] * 18 + [0] * 2})
    assert stats.testable_systems(m, 0.05, 0.95, min_genomes=1) == ["common", "mid", "rare"]
    assert stats.testable_systems(m, 0.05, 0.95, min_genomes=3) == ["mid"]


def test_cooccurrence_counts_and_fisher_value():
    # 2x2 table [[3, 1], [1, 3]]: two-sided Fisher p = 34/70 = 0.4857
    a = [1, 1, 1, 1, 0, 0, 0, 0]
    b = [1, 1, 1, 0, 1, 0, 0, 0]
    cooc = stats.cooccurrence(_matrix({"A": a, "B": b}), 0.05, 0.95)
    row = cooc.iloc[0]
    assert (row.n_both, row.n_a_only, row.n_b_only, row.n_neither) == (3, 1, 1, 3)
    assert row.p_value == pytest.approx(34 / 70)
    assert row.odds_ratio == pytest.approx(9.0)
    assert row.log2_or_haldane == pytest.approx(np.log2(3.5 * 3.5 / (1.5 * 1.5)))
    assert row.association == "ns"


def test_cooccurrence_detects_perfect_exclusion():
    a = [1] * 10 + [0] * 10
    b = [0] * 10 + [1] * 10
    cooc = stats.cooccurrence(_matrix({"A": a, "B": b}), 0.05, 0.95)
    row = cooc.iloc[0]
    assert row.odds_ratio == 0
    assert np.isfinite(row.log2_or_haldane) and row.log2_or_haldane < 0
    assert row.association == "exclusion"


def test_cooccurrence_bh_q_values():
    rng = np.random.default_rng(0)
    m = _matrix({f"S{i}": rng.integers(0, 2, 30).tolist() for i in range(6)})
    cooc = stats.cooccurrence(m, 0.05, 0.95)
    assert len(cooc) == 15  # 6 choose 2
    assert (cooc["q_value"] >= cooc["p_value"]).all()
    ordered = cooc.sort_values("p_value")
    assert ordered["q_value"].is_monotonic_increasing


def test_pair_linkage_distances():
    cooc = pd.DataFrame({"system_a": ["A", "A"], "system_b": ["B", "C"]})
    systems = pd.DataFrame(
        {
            "accession": ["G1", "G1", "G1", "G2", "G2", "G2"],
            "type": ["A", "B", "C", "A", "B", "C"],
            "replicon": ["chr", "chr", "plasmid", "chr", "chr", "chr"],
            "start": [1_000, 3_000, 500, 10_000, 9_000, 900_000],
            "end": [2_000, 4_000, 800, 11_000, 10_500, 901_000],
        }
    )
    out = stats.pair_linkage(cooc, systems, "type", linkage_bp=5_000).set_index("system_b")
    # A-B: gap 1 kb in G1, overlap (0) in G2 -> median 500 bp, linked
    assert out.loc["B", "n_both_same_replicon"] == 2
    assert out.loc["B", "median_gap_bp"] == 500
    assert out.loc["B", "physically_linked"]
    # A-C: different replicons in G1, 889 kb apart in G2 -> not linked
    assert out.loc["C", "n_both_same_replicon"] == 1
    assert not out.loc["C", "physically_linked"]


def test_cooccurrence_with_no_testable_systems():
    cooc = stats.cooccurrence(_matrix({"core": [1, 1, 1]}), 0.05, 0.95)
    assert cooc.empty
    assert "q_value" in cooc.columns
