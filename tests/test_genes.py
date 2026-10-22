import random
from pathlib import Path

import pytest

from pseudodefense import genes


def _random_genome(path: Path, lengths: list[int], seed: int = 1) -> Path:
    # pyrodigal warns below 100 kb of training sequence; random GC-rich DNA is enough to test I/O
    rng = random.Random(seed)
    with path.open("w") as fh:
        for i, n in enumerate(lengths):
            seq = "".join(rng.choices("ACGT", weights=[17, 33, 33, 17], k=n))
            fh.write(f">contig{i + 1} description\n")
            fh.writelines(seq[j : j + 80] + "\n" for j in range(0, n, 80))
    return path


def test_read_fasta_strips_descriptions(tmp_path):
    fasta = tmp_path / "g.fna"
    fasta.write_text(">a desc\nAC\nGT\n>b\nTT\n")
    assert list(genes.read_fasta(fasta)) == [("a", "ACGT"), ("b", "TT")]


def test_predict_genes_outputs_match(tmp_path):
    fasta = _random_genome(tmp_path / "g.fna", [80_000, 30_000])
    faa, gff = tmp_path / "g.faa", tmp_path / "g.gff"
    n = genes.predict_genes(fasta, faa, gff)

    protein_ids = [line[1:].split()[0] for line in faa.read_text().splitlines() if line[:1] == ">"]
    gff_ids = [
        line.split("\t")[8].split(";")[0].removeprefix("ID=")
        for line in gff.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert n == len(protein_ids) > 0
    assert protein_ids == gff_ids
    # IDs are <contig>_<n>, contigs in input order (gembase-like layout for DefenseFinder)
    contigs = [pid.rsplit("_", 1)[0] for pid in protein_ids]
    assert contigs == sorted(contigs)
    assert "*" not in faa.read_text()
    assert gff.read_text().startswith("##gff-version 3\n")


def test_predict_genes_empty_fasta(tmp_path):
    fasta = tmp_path / "empty.fna"
    fasta.write_text("")
    with pytest.raises(ValueError, match="no sequences"):
        genes.predict_genes(fasta, tmp_path / "x.faa", tmp_path / "x.gff")
