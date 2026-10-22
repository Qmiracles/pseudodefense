"""Command-line interface. Snakemake rules call these commands."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from pseudodefense import defense, derep, download, genes, plots, qc, stats, tables

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
    require_ani_ok: Annotated[
        bool, typer.Option(help="Keep only assemblies whose NCBI ANI taxonomy check is OK.")
    ] = True,
) -> None:
    """Select genomes reproducibly from the metadata."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    with summary_jsonl.open() as fh:
        table = download.read_summary(fh)
    if require_ani_ok:
        table = download.filter_taxonomy(table)
    selected = download.select_genomes(table, n, seed)
    out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out, sep="\t", index=False)
    typer.echo(f"Selected {len(selected)} of {len(table)} genomes")


@app.command("download")
def download_cmd(
    selected: Annotated[Path, typer.Option(help="TSV from `select`.")],
    genomes_dir: Annotated[Path, typer.Option(help="Directory for <accession>.fna files.")],
    manifest: Annotated[Path, typer.Option(help="Output manifest.tsv.")],
    batch_size: Annotated[int, typer.Option(help="Genomes per datasets download.")] = 200,
) -> None:
    """Download genome FASTA files (gzip-compressed) and write the manifest."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    table = pd.read_csv(selected, sep="\t")
    download.download_genomes(table, genomes_dir, manifest, batch_size)


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


@app.command("genes")
def genes_cmd(
    fasta: Annotated[Path, typer.Option(help="Genome nucleotide FASTA.")],
    faa: Annotated[Path, typer.Option(help="Output protein FASTA.")],
    gff: Annotated[Path, typer.Option(help="Output GFF3.")],
    translation_table: Annotated[int, typer.Option()] = 11,
) -> None:
    """Predict genes with pyrodigal (single mode)."""
    n = genes.predict_genes(fasta, faa, gff, translation_table)
    typer.echo(f"{fasta.name}: {n} genes")


@app.command("parse-defense")
def parse_defense_cmd(
    accessions: Annotated[Path, typer.Option(help="File with one accession per line.")],
    defensefinder_dir: Annotated[Path, typer.Option(help="Directory with <acc>/ outputs.")],
    genes_dir: Annotated[Path, typer.Option(help="Directory with <acc>.gff files.")],
    systems_out: Annotated[Path, typer.Option(help="Output tidy systems table (TSV).")],
    genes_out: Annotated[Path, typer.Option(help="Output tidy genes table (TSV).")],
) -> None:
    """Combine DefenseFinder outputs into tidy systems and genes tables."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    accs = [line.strip() for line in accessions.read_text().splitlines() if line.strip()]
    systems, genes_table = defense.parse_all(accs, defensefinder_dir, genes_dir)
    systems_out.parent.mkdir(parents=True, exist_ok=True)
    systems.to_csv(systems_out, sep="\t", index=False)
    genes_table.to_csv(genes_out, sep="\t", index=False)
    typer.echo(f"{len(accs)} genomes, {len(systems)} systems, {len(genes_table)} gene rows")


def _read_accessions(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


@app.command("tables")
def tables_cmd(
    systems: Annotated[Path, typer.Option(help="Tidy systems table from `parse-defense`.")],
    accessions: Annotated[Path, typer.Option(help="Genomes in the analysis (one per line).")],
    outdir: Annotated[Path, typer.Option(help="Output directory.")],
) -> None:
    """Presence/absence matrices (type, subtype), counts per genome and subtype table."""
    accs = _read_accessions(accessions)
    table = pd.read_csv(systems, sep="\t")
    outdir.mkdir(parents=True, exist_ok=True)
    for level in tables.LEVELS:
        tables.presence_matrix(table, accs, level).to_csv(outdir / f"matrix_{level}.tsv", sep="\t")
    tables.systems_per_genome(table, accs).to_csv(
        outdir / "systems_per_genome.tsv", sep="\t", index=False
    )
    tables.subtype_table(table, len(accs)).to_csv(outdir / "subtypes.tsv", sep="\t", index=False)


@app.command("stats")
def stats_cmd(
    matrix: Annotated[Path, typer.Option(help="Presence/absence matrix (TSV).")],
    prevalence_out: Annotated[Path, typer.Option(help="Output prevalence table.")],
    cooccurrence_out: Annotated[Path, typer.Option(help="Output co-occurrence table.")],
    min_prevalence: Annotated[float, typer.Option()] = 0.05,
    max_prevalence: Annotated[float, typer.Option()] = 0.95,
    alpha: Annotated[float, typer.Option(help="FDR level (Benjamini-Hochberg).")] = 0.05,
    min_genomes: Annotated[
        int, typer.Option(help="Min genomes with and without a system to test it.")
    ] = 1,
    systems: Annotated[
        Path | None, typer.Option(help="Tidy systems table, to compute pair distances.")
    ] = None,
    level: Annotated[str, typer.Option(help="Matrix level (type or subtype).")] = "type",
    linkage_bp: Annotated[int, typer.Option(help="Max median gap to flag a linked pair.")] = 10_000,
) -> None:
    """Prevalence with Wilson intervals and pairwise Fisher tests with BH correction."""
    m = pd.read_csv(matrix, sep="\t", index_col=0)
    prevalence_out.parent.mkdir(parents=True, exist_ok=True)
    stats.prevalence(m).to_csv(prevalence_out, sep="\t", index=False)
    cooc = stats.cooccurrence(m, min_prevalence, max_prevalence, alpha, min_genomes)
    if systems is not None:
        cooc = stats.pair_linkage(cooc, pd.read_csv(systems, sep="\t"), level, linkage_bp)
    cooc.to_csv(cooccurrence_out, sep="\t", index=False)
    sig = cooc["association"] != "ns"
    msg = f"{len(cooc)} pairs tested, {int(sig.sum())} significant at q < {alpha}"
    if "physically_linked" in cooc:
        msg += f" ({int((sig & cooc['physically_linked']).sum())} physically linked)"
    typer.echo(msg)


@app.command("plots")
def plots_cmd(
    matrix: Annotated[Path, typer.Option(help="Presence/absence matrix (TSV).")],
    prevalence: Annotated[Path, typer.Option(help="Prevalence table.")],
    counts: Annotated[Path, typer.Option(help="systems_per_genome.tsv.")],
    cooccurrence: Annotated[Path, typer.Option(help="Co-occurrence table.")],
    outdir: Annotated[Path, typer.Option(help="Output directory for figures.")],
    taxon: Annotated[str, typer.Option(help="Name used in figure titles.")],
    top_n: Annotated[int, typer.Option()] = 30,
    max_cooc_systems: Annotated[int, typer.Option()] = 40,
    alpha: Annotated[float, typer.Option()] = 0.05,
    cooc_representatives: Annotated[
        Path | None,
        typer.Option(help="representatives_<label>.txt the co-occurrence table was built on."),
    ] = None,
) -> None:
    """Heatmap, prevalence bars, systems-per-genome histogram and co-occurrence matrix."""
    m = pd.read_csv(matrix, sep="\t", index_col=0)
    n = len(m)
    plots.presence_heatmap(
        m, outdir / "presence_heatmap.png", f"Defense systems in {n} {taxon} genomes"
    )
    plots.prevalence_bars(
        pd.read_csv(prevalence, sep="\t"),
        outdir / "prevalence.png",
        f"Most prevalent defense systems ({taxon})",
        top_n,
    )
    plots.systems_histogram(
        pd.read_csv(counts, sep="\t"),
        outdir / "systems_per_genome.png",
        f"Defense systems per genome ({taxon}, n = {n})",
    )
    cooc_title = f"Co-occurrence of defense systems ({n} genomes)"
    if cooc_representatives is not None:
        n_reps = len(_read_accessions(cooc_representatives))
        ani = cooc_representatives.stem.removeprefix("representatives_ani")
        cooc_title = (
            f"Co-occurrence of defense systems ({n_reps} genomes, one per {ani} % ANI cluster)"
        )
    plots.cooccurrence_heatmap(
        pd.read_csv(cooccurrence, sep="\t"),
        outdir / "cooccurrence.png",
        cooc_title,
        alpha,
        max_cooc_systems,
    )


@app.command("dereplicate")
def dereplicate_cmd(
    dist: Annotated[Path, typer.Option(help="All-vs-all `mash dist` output (.tsv or .tsv.gz).")],
    manifest: Annotated[Path, typer.Option(help="manifest.tsv.")],
    qc_report: Annotated[Path, typer.Option(help="qc_report.tsv.")],
    threshold: Annotated[list[float], typer.Option(help="ANI threshold (%); repeatable.")],
    outdir: Annotated[Path, typer.Option(help="Output directory.")],
) -> None:
    """Cluster near-identical genomes and pick one representative per cluster."""
    d = derep.read_mash_dist(dist)
    clusters = derep.cluster_genomes(d, threshold)
    ranks = derep.quality_order(pd.read_csv(manifest, sep="\t"), pd.read_csv(qc_report, sep="\t"))
    clusters = derep.choose_representatives(clusters, ranks, threshold)
    outdir.mkdir(parents=True, exist_ok=True)
    clusters.to_csv(outdir / "clusters.tsv", sep="\t", index=False)
    summ = derep.summary(clusters, threshold)
    summ.to_csv(outdir / "summary.tsv", sep="\t", index=False)
    for t in threshold:
        label = derep.threshold_label(t)
        reps = clusters.loc[clusters[f"rep_{label}"], "accession"]
        (outdir / f"representatives_{label}.txt").write_text("".join(f"{a}\n" for a in reps))
    typer.echo(summ.to_string(index=False))


@app.command("stats-derep")
def stats_derep_cmd(
    matrix: Annotated[Path, typer.Option(help="Presence/absence matrix, all genomes.")],
    systems: Annotated[Path, typer.Option(help="Tidy systems table.")],
    representatives: Annotated[
        list[Path], typer.Option(help="representatives_<label>.txt; repeatable.")
    ],
    outdir: Annotated[Path, typer.Option(help="Output directory.")],
    min_prevalence: Annotated[float, typer.Option()] = 0.05,
    max_prevalence: Annotated[float, typer.Option()] = 0.95,
    alpha: Annotated[float, typer.Option()] = 0.05,
    min_genomes: Annotated[int, typer.Option()] = 1,
    level: Annotated[str, typer.Option()] = "type",
    linkage_bp: Annotated[int, typer.Option()] = 10_000,
) -> None:
    """Prevalence for all genomes and each dereplicated set; co-occurrence per set."""
    m = pd.read_csv(matrix, sep="\t", index_col=0)
    sys_table = pd.read_csv(systems, sep="\t")
    outdir.mkdir(parents=True, exist_ok=True)
    sets = {"all": list(m.index)}
    for path in representatives:
        sets[path.stem.removeprefix("representatives_")] = _read_accessions(path)
    prev = []
    for name, accs in sets.items():
        p = stats.prevalence(m.loc[accs])
        p.insert(0, "genome_set", name)
        prev.append(p)
        if name != "all":
            sub = sys_table[sys_table["accession"].isin(accs)]
            cooc = stats.cooccurrence(
                m.loc[accs], min_prevalence, max_prevalence, alpha, min_genomes
            )
            cooc = stats.pair_linkage(cooc, sub, level, linkage_bp)
            cooc.to_csv(outdir / f"cooccurrence_{name}.tsv", sep="\t", index=False)
            n_sig = int((cooc["association"] != "ns").sum())
            typer.echo(f"{name}: {len(accs)} genomes, {len(cooc)} pairs, {n_sig} with q < {alpha}")
    pd.concat(prev, ignore_index=True).to_csv(
        outdir / "prevalence_by_genome_set.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    app()
