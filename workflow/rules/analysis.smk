rule tables:
    """Presence/absence matrices, systems per genome and subtype table."""
    input:
        systems=rules.parse_defense.output.systems,
        accessions=f"{OUT}/qc/passed_accessions.txt",
    output:
        matrix_type=f"{OUT}/tables/matrix_type.tsv",
        matrix_subtype=f"{OUT}/tables/matrix_subtype.tsv",
        counts=f"{OUT}/tables/systems_per_genome.tsv",
        subtypes=f"{OUT}/tables/subtypes.tsv",
    log:
        f"{OUT}/logs/tables.log",
    shell:
        "pseudodefense tables --systems {input.systems} --accessions {input.accessions} "
        "--outdir {OUT}/tables > {log} 2>&1"


rule stats:
    """Prevalence (Wilson 95 % CI) and co-occurrence (Fisher + Benjamini-Hochberg)."""
    input:
        matrix=f"{OUT}/tables/matrix_{config['analysis']['level']}.tsv",
        systems=rules.parse_defense.output.systems,
    output:
        prevalence=f"{OUT}/stats/prevalence.tsv",
        cooccurrence=f"{OUT}/stats/cooccurrence.tsv",
    params:
        min_prev=config["analysis"]["cooccurrence_min_prevalence"],
        max_prev=config["analysis"]["cooccurrence_max_prevalence"],
        alpha=config["analysis"]["fdr_alpha"],
        min_genomes=config["analysis"]["cooccurrence_min_genomes"],
        level=config["analysis"]["level"],
        linkage_bp=config["analysis"]["linkage_bp"],
    log:
        f"{OUT}/logs/stats.log",
    shell:
        "pseudodefense stats --matrix {input.matrix} --systems {input.systems} "
        "--prevalence-out {output.prevalence} --cooccurrence-out {output.cooccurrence} "
        "--min-prevalence {params.min_prev} --max-prevalence {params.max_prev} "
        "--min-genomes {params.min_genomes} --alpha {params.alpha} "
        "--level {params.level} --linkage-bp {params.linkage_bp} > {log} 2>&1"


rule stats_derep:
    """Prevalence for all genomes and each dereplicated set; co-occurrence per set."""
    input:
        matrix=f"{OUT}/tables/matrix_{config['analysis']['level']}.tsv",
        systems=rules.parse_defense.output.systems,
        reps=rules.dereplicate.output.reps,
    output:
        prevalence=f"{OUT}/stats/prevalence_by_genome_set.tsv",
        cooccurrence=expand(f"{OUT}/stats/cooccurrence_{{label}}.tsv", label=DEREP_LABELS),
    params:
        reps=lambda wc, input: " ".join(f"--representatives {p}" for p in input.reps),
        min_prev=config["analysis"]["cooccurrence_min_prevalence"],
        max_prev=config["analysis"]["cooccurrence_max_prevalence"],
        alpha=config["analysis"]["fdr_alpha"],
        min_genomes=config["analysis"]["cooccurrence_min_genomes"],
        level=config["analysis"]["level"],
        linkage_bp=config["analysis"]["linkage_bp"],
    log:
        f"{OUT}/logs/stats_derep.log",
    shell:
        "pseudodefense stats-derep --matrix {input.matrix} --systems {input.systems} "
        "{params.reps} --outdir {OUT}/stats --min-prevalence {params.min_prev} "
        "--max-prevalence {params.max_prev} --alpha {params.alpha} "
        "--min-genomes {params.min_genomes} --level {params.level} "
        "--linkage-bp {params.linkage_bp} > {log} 2>&1"


COOC_PLOT_SET = config["analysis"]["cooccurrence_plot_set"]
if COOC_PLOT_SET != "all" and COOC_PLOT_SET not in DEREP_LABELS:
    raise ValueError(
        f"analysis.cooccurrence_plot_set = {COOC_PLOT_SET!r}; use 'all' or one of {DEREP_LABELS}"
    )


rule plots:
    """Figures (PNG + SVG). Co-occurrence is drawn for `cooccurrence_plot_set`."""
    input:
        matrix=f"{OUT}/tables/matrix_{config['analysis']['level']}.tsv",
        prevalence=rules.stats.output.prevalence,
        counts=rules.tables.output.counts,
        cooccurrence=(
            rules.stats.output.cooccurrence
            if COOC_PLOT_SET == "all"
            else f"{OUT}/stats/cooccurrence_{COOC_PLOT_SET}.tsv"
        ),
        cooc_reps=(
            [] if COOC_PLOT_SET == "all" else f"{OUT}/derep/representatives_{COOC_PLOT_SET}.txt"
        ),
    output:
        heatmap=f"{OUT}/figures/presence_heatmap.png",
        prevalence=f"{OUT}/figures/prevalence.png",
        histogram=f"{OUT}/figures/systems_per_genome.png",
        cooccurrence=f"{OUT}/figures/cooccurrence.png",
    params:
        taxon=config["taxon"],
        top_n=config["analysis"]["prevalence_top_n"],
        max_cooc=config["analysis"]["cooccurrence_plot_max_systems"],
        alpha=config["analysis"]["fdr_alpha"],
        cooc_reps=lambda wc, input: (
            f"--cooc-representatives {input.cooc_reps}" if input.cooc_reps else ""
        ),
    log:
        f"{OUT}/logs/plots.log",
    shell:
        "pseudodefense plots --matrix {input.matrix} --prevalence {input.prevalence} "
        "--counts {input.counts} --cooccurrence {input.cooccurrence} {params.cooc_reps} "
        "--outdir {OUT}/figures --taxon '{params.taxon}' --top-n {params.top_n} "
        "--max-cooc-systems {params.max_cooc} --alpha {params.alpha} > {log} 2>&1"
