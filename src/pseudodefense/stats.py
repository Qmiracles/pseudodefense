"""Descriptive statistics on the presence/absence matrix.

Results are descriptive associations across genomes, not causal effects.
Genomes are not independent (shared ancestry, near-identical isolates), so
p-values here overstate the evidence; see the README limitations.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats


def wilson_interval(k: int | np.ndarray, n: int, z: float = 1.959964) -> tuple:
    """Wilson score interval for a binomial proportion ``k / n`` (default 95 %)."""
    k = np.asarray(k, dtype=float)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return np.clip(centre - half, 0, 1), np.clip(centre + half, 0, 1)


def prevalence(matrix: pd.DataFrame) -> pd.DataFrame:
    """Fraction of genomes carrying each system, with a 95 % Wilson interval."""
    n = len(matrix)
    k = matrix.sum(axis=0)
    low, high = wilson_interval(k.to_numpy(), n)
    out = pd.DataFrame(
        {
            "system": k.index,
            "n_genomes": k.to_numpy().astype(int),
            "total_genomes": n,
            "prevalence": k.to_numpy() / n,
            "ci95_low": low,
            "ci95_high": high,
        }
    )
    return out.sort_values(["n_genomes", "system"], ascending=[False, True]).reset_index(drop=True)


def testable_systems(
    matrix: pd.DataFrame, min_prevalence: float, max_prevalence: float, min_genomes: int = 1
) -> list[str]:
    """Systems with prevalence in ``[min_prevalence, max_prevalence]``.

    Additionally each system must be present in at least ``min_genomes`` genomes
    and absent from at least ``min_genomes``. Systems seen in (almost) every
    genome or in (almost) none carry no information about co-occurrence: their
    Fisher tests cannot reach significance, yet they would still inflate the
    number of tests in the Benjamini-Hochberg correction.
    """
    n_present = matrix.sum(axis=0)
    prev = n_present / len(matrix)
    keep = (
        (prev >= min_prevalence)
        & (prev <= max_prevalence)
        & (n_present >= min_genomes)
        & (len(matrix) - n_present >= min_genomes)
    )
    return sorted(prev[keep].index)


def pair_linkage(
    cooc: pd.DataFrame, systems: pd.DataFrame, level: str = "type", linkage_bp: int = 10_000
) -> pd.DataFrame:
    """Add the genomic distance between the two systems of each tested pair.

    For every genome carrying both systems on the same replicon, the gap is the
    smallest distance between any copy of A and any copy of B (0 if they
    overlap). A pair is flagged ``physically_linked`` when the median gap is at
    most ``linkage_bp``: such pairs co-occur because they sit in the same locus
    (or DefenseFinder splits one locus into two systems), not necessarily
    because of an independent association.

    Args:
        cooc: output of :func:`cooccurrence`.
        systems: tidy systems table with ``accession``, ``replicon``, ``start``, ``end``.
        level: column of ``systems`` matching the matrix used for ``cooc``.
        linkage_bp: maximum median gap (bp) to call a pair physically linked.

    """
    loci = {
        key: g[["replicon", "start", "end"]].to_numpy()
        for key, g in systems.groupby(["accession", level])
    }
    accessions = systems["accession"].unique()
    n_same, median_gap = [], []
    for a, b in zip(cooc["system_a"], cooc["system_b"], strict=True):
        gaps = []
        for acc in accessions:
            la, lb = loci.get((acc, a)), loci.get((acc, b))
            if la is None or lb is None:
                continue
            same = [
                max(0, max(sa, sb) - min(ea, eb))
                for ra, sa, ea in la
                for rb, sb, eb in lb
                if ra == rb
            ]
            if same:
                gaps.append(min(same))
        n_same.append(len(gaps))
        median_gap.append(float(np.median(gaps)) if gaps else np.nan)
    out = cooc.copy()
    out["n_both_same_replicon"] = n_same
    out["median_gap_bp"] = median_gap
    out["physically_linked"] = out["median_gap_bp"] <= linkage_bp
    return out


def cooccurrence(
    matrix: pd.DataFrame,
    min_prevalence: float = 0.05,
    max_prevalence: float = 0.95,
    alpha: float = 0.05,
    min_genomes: int = 1,
) -> pd.DataFrame:
    """Pairwise co-occurrence / exclusion with Fisher's exact test and BH correction.

    For each pair of testable systems (A, B) the 2x2 table over genomes is
    ``[[both, A only], [B only, neither]]``. The two-sided Fisher p-values are
    corrected with Benjamini-Hochberg across all pairs tested.

    Returns:
        One row per pair with counts, the sample odds ratio (may be 0 or inf),
        a Haldane-Anscombe corrected log2 odds ratio (+0.5 per cell, finite, for
        plotting; with small counts it is biased and can even have the opposite
        sign of the sample odds ratio), ``p_value``, ``q_value`` and ``association``
        (``co-occurrence``, ``exclusion`` or ``ns`` at ``q < alpha``).

    """
    systems = testable_systems(matrix, min_prevalence, max_prevalence, min_genomes)
    rows = []
    for a, b in combinations(systems, 2):
        x, y = matrix[a].to_numpy().astype(bool), matrix[b].to_numpy().astype(bool)
        both, a_only, b_only = int((x & y).sum()), int((x & ~y).sum()), int((~x & y).sum())
        neither = int((~x & ~y).sum())
        res = stats.fisher_exact([[both, a_only], [b_only, neither]], alternative="two-sided")
        log2_or = np.log2(((both + 0.5) * (neither + 0.5)) / ((a_only + 0.5) * (b_only + 0.5)))
        rows.append(
            {
                "system_a": a,
                "system_b": b,
                "n_both": both,
                "n_a_only": a_only,
                "n_b_only": b_only,
                "n_neither": neither,
                "odds_ratio": float(res.statistic),
                "log2_or_haldane": float(log2_or),
                "p_value": float(res.pvalue),
            }
        )
    out = pd.DataFrame(
        rows,
        columns=[
            "system_a",
            "system_b",
            "n_both",
            "n_a_only",
            "n_b_only",
            "n_neither",
            "odds_ratio",
            "log2_or_haldane",
            "p_value",
        ],
    )
    out["q_value"] = stats.false_discovery_control(out["p_value"], method="bh") if len(out) else []
    out["association"] = np.select(
        [(out["q_value"] < alpha) & (out["odds_ratio"] > 1), (out["q_value"] < alpha)],
        ["co-occurrence", "exclusion"],
        default="ns",
    )
    return out.sort_values(["q_value", "system_a", "system_b"]).reset_index(drop=True)
