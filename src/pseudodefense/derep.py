"""Dereplication of near-identical genomes from Mash distances.

Mash distance approximates ``1 - ANI``. Genomes are clustered by average
linkage and the tree is cut at one or more ANI thresholds; each cluster keeps
one representative chosen deterministically by assembly quality.
"""

from __future__ import annotations

import gzip
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def accession_from_path(path: str) -> str:
    """``results/genomes/GCF_000006765.1.fna.gz`` -> ``GCF_000006765.1``."""
    name = Path(path).name
    for suffix in (".gz", ".fna", ".fasta", ".fa"):
        name = name.removesuffix(suffix)
    return name


def read_mash_dist(path: Path) -> pd.DataFrame:
    """Square distance matrix from all-vs-all ``mash dist`` output (optionally gzipped).

    The output has columns reference, query, distance, p-value, shared-hashes.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        long = pd.read_csv(
            fh, sep="\t", header=None, names=["ref", "query", "dist", "p", "shared"],
            usecols=["ref", "query", "dist"],
        )  # fmt: skip
    long["ref"] = long["ref"].map(accession_from_path)
    long["query"] = long["query"].map(accession_from_path)
    square = long.pivot(index="query", columns="ref", values="dist")
    labels = sorted(square.index)
    square = square.reindex(index=labels, columns=labels)
    if square.isna().to_numpy().any():
        raise ValueError(f"{path}: distance matrix is incomplete (missing pairs)")
    values = square.to_numpy()
    # mash dist is symmetric up to rounding; enforce exact symmetry for scipy
    values = (values + values.T) / 2
    np.fill_diagonal(values, 0.0)
    return pd.DataFrame(values, index=labels, columns=labels)


def threshold_label(ani: float) -> str:
    """Column-safe label for an ANI threshold, e.g. 99.5 -> ``ani99.5``."""
    return f"ani{ani:g}"


def cluster_genomes(dist: pd.DataFrame, thresholds: Sequence[float]) -> pd.DataFrame:
    """Cluster genomes at each ANI threshold (average linkage on Mash distance).

    Two genomes end in the same cluster at threshold T when the average-linkage
    merge happens at distance <= 1 - T/100.

    Returns:
        One row per accession, one integer cluster column per threshold.

    """
    out = pd.DataFrame({"accession": dist.index})
    if len(dist) == 1:
        for t in thresholds:
            out[threshold_label(t)] = 1
        return out
    tree = linkage(squareform(dist.to_numpy(), checks=False), method="average")
    for t in thresholds:
        out[threshold_label(t)] = fcluster(tree, t=1 - t / 100, criterion="distance")
    return out


def quality_order(manifest: pd.DataFrame, qc_report: pd.DataFrame) -> pd.DataFrame:
    """Rank genomes for representative choice (best first).

    Order: Complete Genome before other levels; higher CheckM completeness
    minus 5 x contamination (missing CheckM ranks last); fewer sequences;
    accession as a deterministic tie-breaker.
    """
    df = manifest[["accession", "assembly_level"]].merge(
        qc_report[["accession", "n_sequences", "checkm_completeness", "checkm_contamination"]],
        on="accession",
    )
    df["is_complete"] = df["assembly_level"].eq("Complete Genome")
    df["checkm_score"] = (df["checkm_completeness"] - 5 * df["checkm_contamination"]).fillna(
        -np.inf
    )
    df = df.sort_values(
        ["is_complete", "checkm_score", "n_sequences", "accession"],
        ascending=[False, False, True, True],
    )
    df["quality_rank"] = np.arange(len(df))
    return df[["accession", "quality_rank"]]


def choose_representatives(
    clusters: pd.DataFrame, ranks: pd.DataFrame, thresholds: Sequence[float]
) -> pd.DataFrame:
    """Add ``rep_<label>`` boolean columns: the best-ranked genome of each cluster."""
    df = clusters.merge(ranks, on="accession", how="left")
    if df["quality_rank"].isna().any():
        raise ValueError("some clustered genomes are missing from the manifest/QC report")
    for t in thresholds:
        label = threshold_label(t)
        best = df.loc[df.groupby(label)["quality_rank"].idxmin(), "accession"]
        df[f"rep_{label}"] = df["accession"].isin(best)
    return df.drop(columns="quality_rank")


def select_validation_pairs(
    dist: pd.DataFrame,
    thresholds: Sequence[float],
    per_threshold: int,
    window: float,
    seed: int,
) -> pd.DataFrame:
    """Sample genome pairs whose Mash ANI lies near each threshold.

    These are the pairs where a Mash error would move a genome across a cluster
    boundary, so they are the ones worth re-measuring with an alignment-based ANI.

    Args:
        dist: square Mash distance matrix (:func:`read_mash_dist`).
        thresholds: ANI thresholds (%).
        per_threshold: pairs sampled per threshold (fewer if not enough exist).
        window: half-width (ANI points) of the band around each threshold.
        seed: random seed for the sample.

    Returns:
        Columns ``accession_a``, ``accession_b``, ``mash_ani``, ``near_threshold``;
        a pair is listed once, under the first threshold that sampled it.

    """
    labels = np.asarray(dist.index)
    iu, ju = np.triu_indices(len(labels), k=1)
    ani = 100 * (1 - dist.to_numpy()[iu, ju])
    rng = np.random.default_rng(seed)
    rows, seen = [], set()
    for t in thresholds:
        idx = np.flatnonzero(np.abs(ani - t) <= window)
        idx = [i for i in idx if (iu[i], ju[i]) not in seen]
        take = rng.choice(idx, size=min(per_threshold, len(idx)), replace=False) if idx else []
        for i in sorted(take):
            seen.add((iu[i], ju[i]))
            rows.append(
                {
                    "accession_a": labels[iu[i]],
                    "accession_b": labels[ju[i]],
                    "mash_ani": round(float(ani[i]), 4),
                    "near_threshold": t,
                }
            )
    return pd.DataFrame(rows, columns=["accession_a", "accession_b", "mash_ani", "near_threshold"])


def read_fastani(path: Path) -> pd.DataFrame:
    """Parse concatenated FastANI output (query, reference, ANI, mapped, total fragments)."""
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=["accession_a", "accession_b", "fastani"])
    df = pd.read_csv(
        path, sep="\t", header=None, names=["query", "ref", "fastani", "mapped", "total"]
    )
    df["accession_a"] = df["query"].map(accession_from_path)
    df["accession_b"] = df["ref"].map(accession_from_path)
    return df[["accession_a", "accession_b", "fastani"]]


def _pair_key(a: pd.Series, b: pd.Series) -> pd.Series:
    """Order-independent key for a genome pair."""
    return a.where(a < b, b) + "|" + b.where(a < b, a)


def compare_ani(pairs: pd.DataFrame, fastani: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare Mash and FastANI estimates on the validation pairs.

    Returns:
        The per-pair table (``fastani`` is NaN when FastANI reported nothing) and a
        per-threshold summary: pairs, pairs with FastANI, mean and median absolute
        difference (FastANI - Mash) and the fraction of pairs that fall on the same
        side of the threshold with both methods.

    """
    fa = fastani.copy()
    fa["pair"] = _pair_key(fa["accession_a"], fa["accession_b"])
    out = pairs.copy()
    out["pair"] = _pair_key(out["accession_a"], out["accession_b"])
    out = out.merge(fa[["pair", "fastani"]].drop_duplicates("pair"), on="pair", how="left")
    out = out.drop(columns="pair")
    out["difference"] = out["fastani"] - out["mash_ani"]
    rows = []
    for t, g in out.groupby("near_threshold", sort=False):
        ok = g.dropna(subset=["fastani"])
        same = (ok["mash_ani"] >= t) == (ok["fastani"] >= t)
        rows.append(
            {
                "threshold_ani": t,
                "n_pairs": len(g),
                "n_with_fastani": len(ok),
                "mean_difference": round(ok["difference"].mean(), 4) if len(ok) else np.nan,
                "median_abs_difference": round(ok["difference"].abs().median(), 4)
                if len(ok)
                else np.nan,
                "same_side_fraction": round(same.mean(), 3) if len(ok) else np.nan,
            }
        )
    return out, pd.DataFrame(rows)


def summary(clusters: pd.DataFrame, thresholds: Sequence[float]) -> pd.DataFrame:
    """Count clusters, the largest cluster size and singletons per threshold."""
    rows = []
    for t in thresholds:
        sizes = clusters[threshold_label(t)].value_counts()
        rows.append(
            {
                "threshold_ani": t,
                "n_genomes": len(clusters),
                "n_clusters": len(sizes),
                "largest_cluster": int(sizes.max()),
                "singletons": int((sizes == 1).sum()),
            }
        )
    return pd.DataFrame(rows)
