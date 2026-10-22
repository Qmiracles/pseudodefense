import logging
from pathlib import Path

import pandas as pd
import pytest

from pseudodefense import defense

DATA = Path(__file__).parent / "data"
ACC = "GCF_012971705.1"
PREFIX = DATA / "defensefinder" / ACC / f"{ACC}_defense_finder"
GFF = DATA / "genes" / f"{ACC}.gff"


@pytest.fixture
def parsed():
    return defense.parse_genome(
        ACC, Path(f"{PREFIX}_systems.tsv"), Path(f"{PREFIX}_genes.tsv"), GFF
    )


def test_read_gff_coords():
    coords = defense.read_gff_coords(GFF)
    assert len(coords) == 19
    row = coords.set_index("protein_id").loc["NZ_CP052759.1_1051"]
    assert row["contig"] == "NZ_CP052759.1"
    assert row["start"] < row["end"]
    assert row["strand"] in {"+", "-"}


def test_systems_table_shape_and_duplicate_removed(parsed, caplog):
    systems, _ = parsed
    # 15 data lines in the raw file, one of them an exact duplicate
    assert len(systems) == 14
    assert systems["sys_id"].is_unique
    assert list(systems.columns) == defense.SYSTEMS_COLUMNS


def test_duplicate_removal_is_logged(caplog):
    with caplog.at_level(logging.WARNING):
        defense.parse_genome(ACC, Path(f"{PREFIX}_systems.tsv"), Path(f"{PREFIX}_genes.tsv"), GFF)
    assert "1 exact duplicate" in caplog.text


def test_plasmid_system_keeps_its_replicon(parsed):
    systems, _ = parsed
    plasmid = systems[systems["replicon"] == "NZ_CP052760.1"]
    assert plasmid["type"].tolist() == ["SanaTA"]


def test_system_span_covers_its_genes(parsed):
    systems, genes = parsed
    pycsar = systems.set_index("sys_id").loc["NZ_CP052759.1_Pycsar_TIR_6"]
    assert pycsar["n_genes"] == 2
    member = genes[genes["sys_id"] == "NZ_CP052759.1_Pycsar_TIR_6"]
    assert pycsar["start"] == member["start"].min()
    assert pycsar["end"] == member["end"].max()


def test_shared_protein_appears_in_both_systems(parsed):
    _, genes = parsed
    shared = genes[genes["protein_id"] == "NZ_CP052759.1_1237"]
    assert set(shared["sys_id"]) == {
        "NZ_CP052759.1_D12_RM_Type_II_PDDEXK_12",
        "NZ_CP052759.1_DNA_met_RM_Type_II_PDDEXK_13",
    }


def test_genes_table_types(parsed):
    _, genes = parsed
    assert list(genes.columns) == defense.GENES_COLUMNS
    assert len(genes) == 20
    assert genes["i_eval"].dtype == float
    assert genes["start"].dtype == int


def test_missing_column_fails_loudly(tmp_path):
    bad = tmp_path / "systems.tsv"
    pd.read_csv(f"{PREFIX}_systems.tsv", sep="\t").drop(columns="subtype").to_csv(
        bad, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="format may have changed"):
        defense.parse_genome(ACC, bad, Path(f"{PREFIX}_genes.tsv"), GFF)


def test_protein_missing_from_gff_fails(tmp_path):
    gff = tmp_path / "empty.gff"
    gff.write_text("##gff-version 3\n")
    with pytest.raises(ValueError, match="not found in GFF"):
        defense.parse_genome(ACC, Path(f"{PREFIX}_systems.tsv"), Path(f"{PREFIX}_genes.tsv"), gff)


def test_genome_without_systems(tmp_path):
    systems = tmp_path / "s.tsv"
    genes = tmp_path / "g.tsv"
    systems.write_text("\t".join(pd.read_csv(f"{PREFIX}_systems.tsv", sep="\t").columns) + "\n")
    genes.write_text("\t".join(pd.read_csv(f"{PREFIX}_genes.tsv", sep="\t").columns) + "\n")
    s, g = defense.parse_genome(ACC, systems, genes, GFF)
    assert s.empty
    assert g.empty


def test_parse_all_uses_pipeline_layout():
    systems, genes = defense.parse_all([ACC], DATA / "defensefinder", DATA / "genes")
    assert set(systems["accession"]) == {ACC}
    assert len(genes) == 20
