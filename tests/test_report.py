from pathlib import Path

import pandas as pd

from pseudodefense import report


def _tsv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _run_dir(tmp_path: Path) -> Path:
    r = tmp_path / "run"
    (r / "figures").mkdir(parents=True)
    (r / "figures/prevalence.png").write_bytes(b"\x89PNG fake")
    (r / "versions.txt").write_text("run_date: 2026-10-22\npseudodefense: 0.1.0\n")
    _tsv(r / "metadata/manifest.tsv", pd.DataFrame(
        {"accession": ["A", "B", "C"], "assembly_level": ["Complete Genome"] * 3}))  # fmt: skip
    _tsv(r / "qc/qc_report.tsv", pd.DataFrame(
        {"accession": ["A", "B", "C"], "passed": [True, True, False],
         "reasons": ["", "", "length<5500000"]}))  # fmt: skip
    _tsv(r / "tables/defense_systems.tsv", pd.DataFrame(
        {"accession": ["A", "A", "B"], "type": ["Gabija", "Dag", "Dag"],
         "subtype": ["Gabija", "Dag", "Dag"]}))  # fmt: skip
    _tsv(r / "tables/systems_per_genome.tsv", pd.DataFrame(
        {"accession": ["A", "B"], "n_systems": [2, 1]}))  # fmt: skip
    _tsv(r / "stats/prevalence_by_genome_set.tsv", pd.DataFrame(
        {"genome_set": ["all", "all", "ani99.5", "ani99.5"],
         "system": ["Dag", "Gabija", "Dag", "Gabija"], "n_genomes": [2, 1, 1, 1],
         "total_genomes": [2, 2, 1, 1], "prevalence": [1.0, 0.5, 1.0, 1.0]}))  # fmt: skip
    cooc = pd.DataFrame(
        {"system_a": ["Gabija"], "system_b": ["Wadjet"], "n_both": [5], "odds_ratio": [4.0],
         "q_value": [0.01], "association": ["co-occurrence"], "physically_linked": [False],
         "median_gap_bp": [1_500_000.0]})  # fmt: skip
    _tsv(r / "stats/cooccurrence.tsv", cooc)
    _tsv(r / "stats/cooccurrence_ani99.5.tsv", cooc)
    _tsv(r / "derep/summary.tsv", pd.DataFrame(
        {"threshold_ani": [99.5], "n_genomes": [2], "n_clusters": [1],
         "largest_cluster": [2], "singletons": [0]}))  # fmt: skip
    _tsv(r / "derep/clusters.tsv", pd.DataFrame(
        {"accession": ["A", "B"], "ani99.5": [1, 1], "rep_ani99.5": [True, False]}))  # fmt: skip
    return r


def test_md_table_and_format():
    df = pd.DataFrame({"a": [1.0, 0.000131], "b": ["x", "y"]})
    assert report.md_table(df) == "| a | b |\n| --- | --- |\n| 1 | x |\n| 0.00013 | y |"


def test_build_report_markdown_and_html(tmp_path):
    r = _run_dir(tmp_path)
    blocks = report.build_report(r, "Pseudomonas aeruginosa", "ani99.5", {"PAO1": "A", "X": "Z"})
    md = report.to_markdown(blocks)
    assert "3 genomas evaluados, 2 aprobados y 1 descartados" in md
    assert "Dag (100.0 %)" in md  # near-universal type flagged
    assert "cluster de 2 (representante)" in md  # reference genome located
    assert "| X | Z | no |" in md  # reference genome absent from the run
    assert "Gabija | Wadjet" in md  # unlinked significant pair listed
    assert "![Prevalencia de los tipos mas frecuentes](figures/prevalence.png)" in md
    html = report.to_html(blocks, r)
    assert html.startswith("<!doctype html>")
    assert "data:image/png;base64," in html
