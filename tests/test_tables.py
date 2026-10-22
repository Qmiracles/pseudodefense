import pandas as pd
import pytest

from pseudodefense import tables

SYSTEMS = pd.DataFrame(
    {
        "accession": ["G1", "G1", "G1", "G2"],
        "sys_id": ["s1", "s2", "s3", "s4"],
        "type": ["Gabija", "RM_Type_I", "RM_Type_I", "Gabija"],
        "subtype": ["Gabija", "RM_Type_I", "RM_Type_I", "Gabija"],
    }
)
ACCS = ["G1", "G2", "G3"]  # G3 has no detected systems


def test_presence_matrix_keeps_genomes_without_systems():
    m = tables.presence_matrix(SYSTEMS, ACCS)
    assert m.index.tolist() == ACCS
    assert m.columns.tolist() == ["Gabija", "RM_Type_I"]
    assert m.loc["G3"].sum() == 0


def test_presence_matrix_is_binary():
    m = tables.presence_matrix(SYSTEMS, ACCS)
    assert m.loc["G1", "RM_Type_I"] == 1  # two RM_Type_I systems -> presence 1


def test_presence_matrix_rejects_unknown_genomes():
    with pytest.raises(ValueError, match="not in the accession list"):
        tables.presence_matrix(SYSTEMS, ["G1"])


def test_presence_matrix_rejects_bad_level():
    with pytest.raises(ValueError, match="level"):
        tables.presence_matrix(SYSTEMS, ACCS, level="family")


def test_systems_per_genome_counts_and_zeros():
    counts = tables.systems_per_genome(SYSTEMS, ACCS).set_index("accession")
    assert counts.loc["G1", "n_systems"] == 3
    assert counts.loc["G1", "n_types"] == 2
    assert counts.loc["G3"].tolist() == [0, 0, 0]


def test_subtype_table_fraction_uses_all_genomes():
    sub = tables.subtype_table(SYSTEMS, n_genomes=3).set_index("subtype")
    assert sub.loc["Gabija", "n_genomes"] == 2
    assert sub.loc["Gabija", "genome_fraction"] == pytest.approx(2 / 3)
    assert sub.loc["RM_Type_I", "n_systems"] == 2
