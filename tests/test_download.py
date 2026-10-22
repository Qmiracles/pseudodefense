import json
from pathlib import Path

import pandas as pd
import pytest

from pseudodefense import download

# Reduced from a real `datasets summary genome ... --as-json-lines` record (datasets 18.39.0)
REPORT = {
    "accession": "GCF_000750905.1",
    "organism": {
        "infraspecific_names": {"strain": "PSE305"},
        "organism_name": "Pseudomonas aeruginosa",
        "tax_id": 287,
    },
    "assembly_info": {
        "assembly_level": "Chromosome",
        "assembly_name": "PSE305",
        "release_date": "2014-08-28",
    },
    "assembly_stats": {
        "gc_percent": 66,
        "number_of_contigs": 144,
        "total_sequence_length": "6762448",
    },
    "checkm_info": {"completeness": 98.06, "contamination": 1.24},
    "average_nucleotide_identity": {
        "best_ani_match": {"ani": 99.43, "organism_name": "Pseudomonas aeruginosa"},
        "taxonomy_check_status": "OK",
    },
}


def test_parse_report_flattens_fields():
    row = download.parse_report(REPORT)
    assert row["accession"] == "GCF_000750905.1"
    assert row["strain"] == "PSE305"
    assert row["assembly_level"] == "Chromosome"
    assert row["total_length"] == 6762448.0
    assert row["ncbi_contigs"] == 144.0
    assert row["checkm_contamination"] == 1.24
    assert row["ani_taxonomy_check"] == "OK"
    assert row["ani_best_value"] == 99.43


def test_parse_report_tolerates_missing_optional_sections():
    minimal = {k: REPORT[k] for k in download.REQUIRED_SECTIONS}
    row = download.parse_report(minimal)
    assert row["checkm_completeness"] is None
    assert row["ani_taxonomy_check"] is None


def test_parse_report_fails_on_unexpected_format():
    camel_case = {"accession": "GCF_1", "organism": {}, "assemblyInfo": {}, "assemblyStats": {}}
    with pytest.raises(ValueError, match="format may have changed"):
        download.parse_report(camel_case)


def test_read_summary_drops_duplicates_and_blank_lines():
    lines = [json.dumps(REPORT), "", json.dumps(REPORT)]
    table = download.read_summary(lines)
    assert len(table) == 1


def test_filter_taxonomy_keeps_ok_only():
    table = pd.DataFrame(
        {"accession": ["A", "B", "C"], "ani_taxonomy_check": ["OK", "Inconclusive", None]}
    )
    assert download.filter_taxonomy(table)["accession"].tolist() == ["A"]


def _table(n):
    return pd.DataFrame({"accession": [f"GCF_{i:09d}.1" for i in range(n)]})


def test_select_genomes_is_reproducible_and_order_independent():
    table = _table(50)
    a = download.select_genomes(table, 10, seed=1)
    b = download.select_genomes(table.sample(frac=1, random_state=7), 10, seed=1)
    assert len(a) == 10
    assert a["accession"].tolist() == b["accession"].tolist()


def test_select_genomes_all_when_n_is_none_or_too_large():
    table = _table(5)
    assert len(download.select_genomes(table, None, seed=1)) == 5
    assert len(download.select_genomes(table, 99, seed=1)) == 5


def test_summary_command_flags():
    cmd = download.summary_command(
        "Pseudomonas aeruginosa", ["complete", "chromosome"], "RefSeq", True
    )
    assert cmd[:5] == ["datasets", "summary", "genome", "taxon", "Pseudomonas aeruginosa"]
    assert "complete,chromosome" in cmd
    assert "--exclude-atypical" in cmd


def test_collect_fastas_finds_one_file_per_accession(tmp_path: Path):
    for acc in ["GCF_1.1", "GCF_2.1"]:
        d = tmp_path / acc
        d.mkdir()
        (d / f"{acc}_ASM_genomic.fna").write_text(">x\nACGT\n")
    found = download.collect_fastas(tmp_path, ["GCF_1.1", "GCF_2.1"])
    assert set(found) == {"GCF_1.1", "GCF_2.1"}


def test_collect_fastas_fails_loudly_on_missing_accession(tmp_path: Path):
    with pytest.raises(RuntimeError, match="GCF_9.1"):
        download.collect_fastas(tmp_path, ["GCF_9.1"])
