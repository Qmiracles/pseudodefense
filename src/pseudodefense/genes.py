"""Gene prediction with pyrodigal (Prodigal single mode).

Writes a protein FASTA (for DefenseFinder) and a GFF3 (needed by PADLOC) per
genome. Protein IDs follow Prodigal's convention ``<contig>_<n>`` and are
written in genomic order, contig by contig.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pyrodigal


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(id, sequence)`` pairs from a nucleotide FASTA file."""
    seq_id: str | None = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if seq_id is not None:
        yield seq_id, "".join(chunks)


def predict_genes(fasta: Path, faa: Path, gff: Path, translation_table: int = 11) -> int:
    """Predict genes in one genome and write proteins (``faa``) and GFF3 (``gff``).

    The gene finder is trained once on all sequences of the genome together
    (chromosome and plasmids), as Prodigal does in single mode.

    Returns:
        Number of predicted genes.

    """
    records = list(read_fasta(fasta))
    if not records:
        raise ValueError(f"{fasta}: no sequences found")
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*(seq.encode() for _, seq in records), translation_table=translation_table)

    faa.parent.mkdir(parents=True, exist_ok=True)
    gff.parent.mkdir(parents=True, exist_ok=True)
    n_genes = 0
    with faa.open("w") as faa_fh, gff.open("w") as gff_fh:
        gff_fh.write("##gff-version 3\n")
        for seq_id, seq in records:
            genes = finder.find_genes(seq.encode())
            genes.write_translations(faa_fh, sequence_id=seq_id, include_stop=False)
            genes.write_gff(gff_fh, sequence_id=seq_id, header=False)
            n_genes += len(genes)
    return n_genes
