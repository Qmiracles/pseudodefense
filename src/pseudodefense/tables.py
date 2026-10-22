"""Summary tables built from the tidy DefenseFinder systems table.

All tables are indexed by the full list of genomes that passed QC, so genomes
without any detected system are kept (as zeros), not silently dropped.
"""

from __future__ import annotations

import pandas as pd

LEVELS = ("type", "subtype")


def presence_matrix(
    systems: pd.DataFrame, accessions: list[str], level: str = "type"
) -> pd.DataFrame:
    """Genome x system presence/absence matrix (0/1).

    Args:
        systems: tidy systems table (``defense.SYSTEMS_COLUMNS``).
        accessions: every genome in the analysis, including those with no systems.
        level: ``"type"`` or ``"subtype"`` (DefenseFinder classification).

    Returns:
        Integer DataFrame, rows = accessions (in the given order), columns =
        systems sorted alphabetically.

    """
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}, got {level!r}")
    unknown = set(systems["accession"]) - set(accessions)
    if unknown:
        raise ValueError(
            f"systems table has genomes not in the accession list: {sorted(unknown)[:5]}"
        )
    matrix = (
        pd.crosstab(systems["accession"], systems[level])
        .clip(upper=1)
        .reindex(index=accessions, fill_value=0)
    )
    matrix = matrix.reindex(columns=sorted(matrix.columns))
    matrix.index.name = "accession"
    matrix.columns.name = None
    return matrix.astype(int)


def systems_per_genome(systems: pd.DataFrame, accessions: list[str]) -> pd.DataFrame:
    """Per-genome counts: systems (DefenseFinder rows), distinct types and subtypes."""
    grouped = systems.groupby("accession")
    counts = pd.DataFrame(
        {
            "n_systems": grouped.size(),
            "n_types": grouped["type"].nunique(),
            "n_subtypes": grouped["subtype"].nunique(),
        }
    )
    counts = counts.reindex(accessions, fill_value=0).astype(int)
    counts.index.name = "accession"
    return counts.reset_index()


def subtype_table(systems: pd.DataFrame, n_genomes: int) -> pd.DataFrame:
    """Each type/subtype with the number of systems and of genomes carrying it."""
    out = (
        systems.groupby(["type", "subtype"])
        .agg(n_systems=("sys_id", "size"), n_genomes=("accession", "nunique"))
        .reset_index()
    )
    out["genome_fraction"] = out["n_genomes"] / n_genomes
    return out.sort_values(
        ["n_genomes", "type", "subtype"], ascending=[False, True, True]
    ).reset_index(drop=True)
