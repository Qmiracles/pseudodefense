"""Command-line interface. Snakemake rules call these commands."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from pseudodefense import download, qc

app = typer.Typer(help="Survey of anti-phage defense systems in Pseudomonas genomes.")


@app.command()
def summary(
    taxon: Annotated[str, typer.Option(help="NCBI taxon name or taxid.")],
    out: Annotated[Path, typer.Option(help="Output JSON Lines file.")],
    assembly_level: Annotated[list[str], typer.Option(help="Repeat for several levels.")],
    assembly_source: Annotated[str, typer.Option()] = "RefSeq",
    exclude_atypical: Annotated[bool, typer.Option()] = True,
) -> None:
    """Fetch NCBI assembly metadata for a taxon."""
    download.fetch_summary(taxon, assembly_level, assembly_source, exclude_atypical, out)


@app.command()
def select(
    summary_jsonl: Annotated[Path, typer.Option("--summary", help="JSON Lines from `summary`.")],
    out: Annotated[Path, typer.Option(help="Output TSV with selected genomes.")],
    n: Annotated[int | None, typer.Option(help="Number of genomes; omit for all.")] = None,
    seed: Annotated[int, typer.Option(help="Random seed for sampling.")] = 42,
) -> None:
    """Select genomes reproducibly from the metadata."""
    with summary_jsonl.open() as fh:
        table = download.read_summary(fh)
    selected = download.select_genomes(table, n, seed)
    out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out, sep="\t", index=False)
    typer.echo(f"Selected {len(selected)} of {len(table)} genomes")


@app.command("download")
def download_cmd(
    selected: Annotated[Path, typer.Option(help="TSV from `select`.")],
    genomes_dir: Annotated[Path, typer.Option(help="Directory for <accession>.fna files.")],
    manifest: Annotated[Path, typer.Option(help="Output manifest.tsv.")],
) -> None:
    """Download genome FASTA files and write the manifest."""
    table = pd.read_csv(selected, sep="\t")
    download.download_genomes(table, genomes_dir, manifest)


@app.command("qc")
def qc_cmd(
    manifest: Annotated[Path, typer.Option(help="manifest.tsv from `download`.")],
    report: Annotated[Path, typer.Option(help="Output QC report (TSV).")],
    passed: Annotated[Path, typer.Option(help="Output list of accessions that pass.")],
    min_length: Annotated[int, typer.Option()],
    max_length: Annotated[int, typer.Option()],
    max_sequences: Annotated[int, typer.Option()],
    min_gc: Annotated[float, typer.Option()],
    max_gc: Annotated[float, typer.Option()],
    min_checkm_completeness: Annotated[float, typer.Option()],
    max_checkm_contamination: Annotated[float, typer.Option()],
) -> None:
    """Filter assemblies by size, contigs, GC and CheckM; report what was discarded."""
    thresholds = qc.QCThresholds(
        min_length=min_length,
        max_length=max_length,
        max_sequences=max_sequences,
        min_gc=min_gc,
        max_gc=max_gc,
        min_checkm_completeness=min_checkm_completeness,
        max_checkm_contamination=max_checkm_contamination,
    )
    result = qc.run_qc(pd.read_csv(manifest, sep="\t"), thresholds)
    report.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(report, sep="\t", index=False)
    passed.write_text("".join(f"{a}\n" for a in result.loc[result["passed"], "accession"]))
    typer.echo(qc.qc_summary(result), nl=False)


if __name__ == "__main__":
    app()
