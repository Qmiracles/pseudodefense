rule report:
    """Run report (Markdown + self-contained HTML) from the pipeline outputs."""
    input:
        versions=f"{OUT}/versions.txt",
        manifest=rules.download_genomes.output.manifest,
        qc=f"{OUT}/qc/qc_report.tsv",
        systems=rules.parse_defense.output.systems,
        counts=rules.tables.output.counts,
        prevalence=rules.stats_derep.output.prevalence,
        cooccurrence=rules.stats.output.cooccurrence,
        cooccurrence_derep=rules.stats_derep.output.cooccurrence,
        derep=rules.dereplicate.output.summary,
        clusters=rules.dereplicate.output.clusters,
        ani=rules.ani_validation.output.summary,
        figures=rules.plots.output,
    output:
        md=f"{OUT}/report.md",
        html=f"{OUT}/report.html",
    params:
        taxon=config["taxon"],
        cooc_set=config["analysis"]["cooccurrence_plot_set"],
        refs=" ".join(
            f"--reference {name}={acc}"
            for name, acc in config["report"]["reference_genomes"].items()
        ),
        alpha=config["analysis"]["fdr_alpha"],
        linkage_bp=config["analysis"]["linkage_bp"],
    log:
        f"{OUT}/logs/report.log",
    shell:
        "pseudodefense report --results {OUT} --md-out {output.md} --html-out {output.html} "
        "--taxon '{params.taxon}' --cooc-set {params.cooc_set} {params.refs} "
        "--alpha {params.alpha} --linkage-bp {params.linkage_bp} > {log} 2>&1"
