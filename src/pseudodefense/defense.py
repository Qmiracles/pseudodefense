"""Parse DefenseFinder output into tidy tables.

DefenseFinder (3.0.0) writes, per genome, ``<name>_defense_finder_systems.tsv``
(one row per system) and ``<name>_defense_finder_genes.tsv`` (one row per gene
in a system). Here they are combined with the pyrodigal GFF3 to add nucleotide
coordinates, producing:

* systems table: one row per (genome, system)
* genes table: one row per (genome, system, protein)
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

SYSTEMS_REQUIRED = ["sys_id", "type", "subtype", "activity", "protein_in_syst", "genes_count"]
GENES_REQUIRED = [
    "replicon",
    "hit_id",
    "gene_name",
    "hit_pos",
    "sys_id",
    "hit_status",
    "hit_i_eval",
    "hit_score",
    "hit_profile_cov",
    "hit_seq_cov",
    "type",
    "subtype",
    "activity",
]

SYSTEMS_COLUMNS = [
    "accession",
    "sys_id",
    "replicon",
    "type",
    "subtype",
    "activity",
    "n_genes",
    "proteins",
    "gene_names",
    "start",
    "end",
    "strands",
]
GENES_COLUMNS = [
    "accession",
    "sys_id",
    "replicon",
    "type",
    "subtype",
    "activity",
    "protein_id",
    "gene_name",
    "hit_status",
    "start",
    "end",
    "strand",
    "i_eval",
    "score",
    "profile_cov",
    "seq_cov",
]


def _read_table(path: Path, required: list[str]) -> pd.DataFrame:
    """Read a DefenseFinder TSV and check that the expected columns exist."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: missing columns {missing}; the DefenseFinder output format may have changed"
        )
    return df


def drop_exact_duplicates(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """Remove fully identical rows, logging how many were dropped.

    DefenseFinder 3.0.0 can write the same system line twice in the systems
    table; identical rows carry no extra information.
    """
    deduped = df.drop_duplicates()
    n = len(df) - len(deduped)
    if n:
        log.warning("%s: dropped %d exact duplicate row(s)", label, n)
    return deduped.reset_index(drop=True)


def read_gff_coords(gff: Path) -> pd.DataFrame:
    """Protein coordinates from a pyrodigal GFF3: protein_id, contig, start, end, strand."""
    rows = []
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in fields[8].strip(";").split(";") if "=" in kv)
            rows.append(
                {
                    "protein_id": attrs["ID"],
                    "contig": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "contig", "start", "end", "strand"])


def tidy_genes(accession: str, genes: pd.DataFrame, coords: pd.DataFrame) -> pd.DataFrame:
    """One row per protein per system, with nucleotide coordinates."""
    merged = genes.merge(coords, left_on="hit_id", right_on="protein_id", how="left")
    unmatched = merged.loc[merged["start"].isna(), "hit_id"].tolist()
    if unmatched:
        raise ValueError(f"{accession}: proteins not found in GFF: {unmatched[:5]}")
    out = merged.rename(
        columns={
            "hit_i_eval": "i_eval",
            "hit_score": "score",
            "hit_profile_cov": "profile_cov",
            "hit_seq_cov": "seq_cov",
        }
    )
    out.insert(0, "accession", accession)
    out = out[GENES_COLUMNS].copy()
    for col in ["start", "end"]:
        out[col] = out[col].astype(int)
    for col in ["i_eval", "score", "profile_cov", "seq_cov"]:
        out[col] = out[col].astype(float)
    return out.sort_values(["replicon", "start", "sys_id"]).reset_index(drop=True)


def tidy_systems(accession: str, systems: pd.DataFrame, genes_tidy: pd.DataFrame) -> pd.DataFrame:
    """One row per system, with replicon and the span of its genes.

    ``start``/``end`` are the minimum and maximum coordinates of the system's
    genes; a system crossing the origin of a circular replicon would get a
    span covering the whole replicon (not handled).
    """
    agg = (
        genes_tidy.sort_values("start")
        .groupby("sys_id")
        .agg(
            replicon=("replicon", "first"),
            n_replicons=("replicon", "nunique"),
            proteins=("protein_id", ",".join),
            gene_names=("gene_name", ",".join),
            start=("start", "min"),
            end=("end", "max"),
            strands=("strand", lambda s: ",".join(sorted(set(s)))),
            n_gene_rows=("protein_id", "size"),
        )
        .reset_index()
    )
    if (agg["n_replicons"] > 1).any():
        raise ValueError(f"{accession}: a system spans several replicons")
    merged = systems.merge(agg, on="sys_id", how="left")
    missing = merged.loc[merged["replicon"].isna(), "sys_id"].tolist()
    if missing:
        raise ValueError(f"{accession}: systems without genes in genes table: {missing}")
    merged["n_genes"] = merged["genes_count"].astype(int)
    mismatch = merged.loc[merged["n_genes"] != merged["n_gene_rows"], "sys_id"].tolist()
    if mismatch:
        raise ValueError(f"{accession}: genes_count differs from genes table for {mismatch}")
    merged.insert(0, "accession", accession)
    out = merged[SYSTEMS_COLUMNS].copy()
    out[["start", "end"]] = out[["start", "end"]].astype(int)
    return out.sort_values(["replicon", "start", "sys_id"]).reset_index(drop=True)


def parse_genome(
    accession: str, systems_tsv: Path, genes_tsv: Path, gff: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Tidy systems and genes tables for one genome."""
    systems = drop_exact_duplicates(
        _read_table(systems_tsv, SYSTEMS_REQUIRED), f"{accession} systems"
    )
    genes = drop_exact_duplicates(_read_table(genes_tsv, GENES_REQUIRED), f"{accession} genes")
    if systems["sys_id"].duplicated().any():
        dup = systems.loc[systems["sys_id"].duplicated(), "sys_id"].tolist()
        raise ValueError(f"{accession}: conflicting rows for systems {dup}")
    genes_tidy = tidy_genes(accession, genes, read_gff_coords(gff))
    return tidy_systems(accession, systems, genes_tidy), genes_tidy


def parse_all(
    accessions: list[str], defensefinder_dir: Path, genes_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Tidy tables for many genomes, using the pipeline's directory layout."""
    systems, genes = [], []
    for acc in accessions:
        prefix = defensefinder_dir / acc / f"{acc}_defense_finder"
        s, g = parse_genome(
            acc,
            Path(f"{prefix}_systems.tsv"),
            Path(f"{prefix}_genes.tsv"),
            genes_dir / f"{acc}.gff",
        )
        systems.append(s)
        genes.append(g)
    return (
        pd.concat(systems, ignore_index=True) if systems else pd.DataFrame(columns=SYSTEMS_COLUMNS),
        pd.concat(genes, ignore_index=True) if genes else pd.DataFrame(columns=GENES_COLUMNS),
    )
