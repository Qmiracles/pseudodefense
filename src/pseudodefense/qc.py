"""Assembly quality control.

Statistics are recomputed from the downloaded FASTA (not only taken from NCBI
metadata). CheckM values come from the NCBI assembly report when available.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class QCThresholds:
    """Configurable QC limits (see ``config/config.yaml``)."""

    min_length: int
    max_length: int
    max_sequences: int
    min_gc: float
    max_gc: float
    min_checkm_completeness: float
    max_checkm_contamination: float


@dataclass(frozen=True)
class FastaStats:
    """Basic assembly statistics computed from a FASTA file."""

    length: int
    n_sequences: int
    gc_percent: float


def fasta_stats(path: Path) -> FastaStats:
    """Compute total length, number of sequences and GC% of a nucleotide FASTA.

    GC% is computed over unambiguous bases (A, C, G, T) only.
    """
    length = n_sequences = gc = acgt = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n_sequences += 1
                continue
            seq = line.strip().upper()
            length += len(seq)
            g_c = seq.count("G") + seq.count("C")
            gc += g_c
            acgt += g_c + seq.count("A") + seq.count("T")
    gc_percent = 100 * gc / acgt if acgt else 0.0
    return FastaStats(length=length, n_sequences=n_sequences, gc_percent=gc_percent)


def failure_reasons(
    stats: FastaStats,
    completeness: float | None,
    contamination: float | None,
    t: QCThresholds,
) -> list[str]:
    """Return the list of QC criteria a genome fails (empty list = pass)."""
    reasons = []
    if stats.length < t.min_length:
        reasons.append(f"length<{t.min_length}")
    if stats.length > t.max_length:
        reasons.append(f"length>{t.max_length}")
    if stats.n_sequences > t.max_sequences:
        reasons.append(f"sequences>{t.max_sequences}")
    if stats.gc_percent < t.min_gc:
        reasons.append(f"gc<{t.min_gc}")
    if stats.gc_percent > t.max_gc:
        reasons.append(f"gc>{t.max_gc}")
    if completeness is not None and completeness < t.min_checkm_completeness:
        reasons.append(f"checkm_completeness<{t.min_checkm_completeness}")
    if contamination is not None and contamination > t.max_checkm_contamination:
        reasons.append(f"checkm_contamination>{t.max_checkm_contamination}")
    return reasons


def _optional(value: object) -> float | None:
    return None if pd.isna(value) else float(value)


def run_qc(manifest: pd.DataFrame, thresholds: QCThresholds) -> pd.DataFrame:
    """Evaluate every genome in the manifest and return a QC report.

    The report has one row per genome with the observed statistics, a boolean
    ``passed`` column and a ``reasons`` column (``;``-separated).
    """
    rows = []
    for rec in manifest.itertuples(index=False):
        stats = fasta_stats(Path(rec.fasta))
        completeness = _optional(rec.checkm_completeness)
        contamination = _optional(rec.checkm_contamination)
        reasons = failure_reasons(stats, completeness, contamination, thresholds)
        rows.append(
            {
                "accession": rec.accession,
                "length": stats.length,
                "n_sequences": stats.n_sequences,
                "gc_percent": round(stats.gc_percent, 3),
                "checkm_completeness": completeness,
                "checkm_contamination": contamination,
                "passed": not reasons,
                "reasons": ";".join(reasons),
            }
        )
    return pd.DataFrame(rows)


def qc_summary(report: pd.DataFrame) -> str:
    """Human-readable summary: how many genomes were discarded and why."""
    n_fail = int((~report["passed"]).sum())
    lines = [f"Genomes evaluated: {len(report)}", f"Passed: {len(report) - n_fail}"]
    lines.append(f"Discarded: {n_fail}")
    if n_fail:
        counts = report.loc[~report["passed"], "reasons"].str.split(";").explode().value_counts()
        lines += [f"  {reason}: {count}" for reason, count in counts.items()]
    return "\n".join(lines) + "\n"
