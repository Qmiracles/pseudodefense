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
