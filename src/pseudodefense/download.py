"""Genome selection and download from NCBI with the `datasets` CLI.

Flow: ``datasets summary`` (metadata as JSON Lines) -> ``select_genomes``
(reproducible sampling) -> ``datasets download`` (FASTA) -> ``manifest.tsv``.
"""

from __future__ import annotations

import datetime as dt
import json
import shutil
import subprocess
import tempfile
import zipfile
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

MANIFEST_COLUMNS = [
    "accession",
    "organism",
    "strain",
    "assembly_level",
    "assembly_name",
    "release_date",
    "total_length",
    "ncbi_contigs",
    "gc_percent",
    "checkm_completeness",
    "checkm_contamination",
    "ani_taxonomy_check",
    "ani_best_match",
    "ani_best_value",
    "fasta",
    "download_date",
    "datasets_version",
]


def datasets_version(executable: str = "datasets") -> str:
    """Return the version string reported by the NCBI ``datasets`` CLI."""
    out = subprocess.run(
        [executable, "--version"], check=True, capture_output=True, text=True
    ).stdout
    # Output looks like "datasets version: 18.39.0"
    return out.strip().split()[-1]


def summary_command(
    taxon: str,
    assembly_levels: Sequence[str],
    assembly_source: str,
    exclude_atypical: bool,
) -> list[str]:
    """Build the ``datasets summary genome taxon`` command."""
    cmd = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        taxon,
        "--assembly-level",
        ",".join(assembly_levels),
        "--assembly-source",
        assembly_source,
        "--as-json-lines",
    ]
    if exclude_atypical:
        cmd.append("--exclude-atypical")
    return cmd


def fetch_summary(
    taxon: str,
    assembly_levels: Sequence[str],
    assembly_source: str,
    exclude_atypical: bool,
    out_jsonl: Path,
) -> None:
    """Write NCBI assembly reports for ``taxon`` to ``out_jsonl`` (one JSON per line)."""
    cmd = summary_command(taxon, assembly_levels, assembly_source, exclude_atypical)
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with out_jsonl.open("w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)


def _to_float(value: object) -> float | None:
    return None if value in (None, "") else float(value)


# Sections every report must have. datasets v18 uses snake_case keys
# (assembly_info, assembly_stats, ...); a missing section means the output
# format changed, so parsing stops instead of producing an empty manifest.
REQUIRED_SECTIONS = ("accession", "organism", "assembly_info", "assembly_stats")


def parse_report(record: dict) -> dict:
    """Flatten one NCBI assembly report into a manifest row (without download fields).

    ``checkm_info`` and ``average_nucleotide_identity`` are optional: NCBI does
    not compute them for every assembly.
    """
    missing = [key for key in REQUIRED_SECTIONS if key not in record]
    if missing:
        raise ValueError(
            f"{record.get('accession', '?')}: missing {missing} in NCBI assembly report; "
            "the `datasets` output format may have changed"
        )
    organism = record["organism"]
    info = record["assembly_info"]
    stats = record["assembly_stats"]
    checkm = record.get("checkm_info", {})
    ani = record.get("average_nucleotide_identity", {})
    return {
        "accession": record["accession"],
        "organism": organism.get("organism_name"),
        "strain": organism.get("infraspecific_names", {}).get("strain"),
        "assembly_level": info.get("assembly_level"),
        "assembly_name": info.get("assembly_name"),
        "release_date": info.get("release_date"),
        "total_length": _to_float(stats.get("total_sequence_length")),
        "ncbi_contigs": _to_float(stats.get("number_of_contigs")),
        "gc_percent": _to_float(stats.get("gc_percent")),
        "checkm_completeness": _to_float(checkm.get("completeness")),
        "checkm_contamination": _to_float(checkm.get("contamination")),
        "ani_taxonomy_check": ani.get("taxonomy_check_status"),
        "ani_best_match": ani.get("best_ani_match", {}).get("organism_name"),
        "ani_best_value": _to_float(ani.get("best_ani_match", {}).get("ani")),
    }


def read_summary(lines: Iterable[str]) -> pd.DataFrame:
    """Parse JSON Lines assembly reports into a table, one row per accession."""
    rows = [parse_report(json.loads(line)) for line in lines if line.strip()]
    df = pd.DataFrame(rows, columns=MANIFEST_COLUMNS[:14])
    return df.drop_duplicates("accession").reset_index(drop=True)


def select_genomes(summary: pd.DataFrame, n: int | None, seed: int) -> pd.DataFrame:
    """Pick ``n`` genomes reproducibly (all of them if ``n`` is None).

    The table is sorted by accession before sampling so the result does not
    depend on the order in which NCBI returns the reports.
    """
    ordered = summary.sort_values("accession").reset_index(drop=True)
    if n is None or n >= len(ordered):
        return ordered
    return ordered.sample(n=n, random_state=seed).sort_values("accession").reset_index(drop=True)


def download_genomes(
    selected: pd.DataFrame, genomes_dir: Path, manifest_path: Path
) -> pd.DataFrame:
    """Download FASTA files for ``selected`` accessions and write the manifest.

    Each genome ends up as ``genomes_dir/<accession>.fna``.
    """
    genomes_dir.mkdir(parents=True, exist_ok=True)
    version = datasets_version()
    today = dt.date.today().isoformat()
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("\n".join(selected["accession"]) + "\n")
        zip_path = tmp_path / "ncbi_dataset.zip"
        subprocess.run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                "--inputfile",
                str(acc_file),
                "--include",
                "genome",
                "--filename",
                str(zip_path),
            ],
            check=True,
        )
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp_path)
        fastas = collect_fastas(tmp_path / "ncbi_dataset" / "data", selected["accession"])
        for accession, src in fastas.items():
            shutil.move(src, genomes_dir / f"{accession}.fna")

    manifest = selected.copy()
    manifest["fasta"] = [str(genomes_dir / f"{a}.fna") for a in manifest["accession"]]
    manifest["download_date"] = today
    manifest["datasets_version"] = version
    manifest = manifest[MANIFEST_COLUMNS]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, sep="\t", index=False)
    return manifest


def collect_fastas(data_dir: Path, accessions: Iterable[str]) -> dict[str, Path]:
    """Map each accession to its single genomic FASTA inside an extracted package.

    Raises an error if any accession is missing or has more than one FASTA,
    instead of silently skipping it.
    """
    found: dict[str, Path] = {}
    problems: list[str] = []
    for accession in accessions:
        hits = sorted((data_dir / accession).glob("*_genomic.fna"))
        if len(hits) != 1:
            problems.append(f"{accession}: {len(hits)} FASTA files")
        else:
            found[accession] = hits[0]
    if problems:
        raise RuntimeError("Unexpected NCBI Datasets package layout:\n" + "\n".join(problems))
    return found
