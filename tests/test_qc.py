from pathlib import Path

import pandas as pd
import pytest

from pseudodefense import qc

THRESHOLDS = qc.QCThresholds(
    min_length=8,
    max_length=100,
    max_sequences=2,
    min_gc=40.0,
    max_gc=60.0,
    min_checkm_completeness=95.0,
    max_checkm_contamination=5.0,
)


def _fasta(path: Path, *seqs: str) -> Path:
    path.write_text("".join(f">c{i}\n{s}\n" for i, s in enumerate(seqs)))
    return path


def test_fasta_stats(tmp_path: Path):
    stats = qc.fasta_stats(_fasta(tmp_path / "g.fna", "ACGT", "GGCCNN"))
    assert stats.length == 10
    assert stats.n_sequences == 2
    # 6 G/C over 8 unambiguous bases; N is ignored
    assert stats.gc_percent == pytest.approx(75.0)


def test_fasta_stats_reads_gzip(tmp_path: Path):
    import gzip

    path = tmp_path / "g.fna.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">c0\nACGT\n>c1\nGGCCNN\n")
    stats = qc.fasta_stats(path)
    assert (stats.length, stats.n_sequences) == (10, 2)


def test_failure_reasons_pass():
    stats = qc.FastaStats(length=50, n_sequences=1, gc_percent=50.0)
    assert qc.failure_reasons(stats, 99.0, 1.0, THRESHOLDS) == []


def test_failure_reasons_lists_every_failed_criterion():
    stats = qc.FastaStats(length=5, n_sequences=3, gc_percent=70.0)
    reasons = qc.failure_reasons(stats, 90.0, 10.0, THRESHOLDS)
    assert reasons == [
        "length<8",
        "sequences>2",
        "gc>60.0",
        "checkm_completeness<95.0",
        "checkm_contamination>5.0",
    ]


def test_missing_checkm_is_not_a_failure():
    stats = qc.FastaStats(length=50, n_sequences=1, gc_percent=50.0)
    assert qc.failure_reasons(stats, None, None, THRESHOLDS) == []


def test_run_qc_and_summary(tmp_path: Path):
    good = _fasta(tmp_path / "good.fna", "ACGTACGTGC")
    bad = _fasta(tmp_path / "bad.fna", "GGGG")
    manifest = pd.DataFrame(
        {
            "accession": ["A", "B"],
            "fasta": [str(good), str(bad)],
            "checkm_completeness": [99.0, float("nan")],
            "checkm_contamination": [0.5, float("nan")],
        }
    )
    report = qc.run_qc(manifest, THRESHOLDS)
    assert report["passed"].tolist() == [True, False]
    assert report.loc[1, "reasons"] == "length<8;gc>60.0"
    text = qc.qc_summary(report)
    assert "Discarded: 1" in text
    assert "length<8: 1" in text
