rule genome_summary:
    """Metadata for every assembly of the taxon that passes the NCBI filters."""
    output:
        f"{OUT}/metadata/genome_summary.jsonl",
    params:
        taxon=config["taxon"],
        levels=" ".join(f"--assembly-level {lvl}" for lvl in config["assembly_levels"]),
        source=config["assembly_source"],
        atypical="--exclude-atypical" if config["exclude_atypical"] else "--no-exclude-atypical",
    log:
        f"{OUT}/logs/genome_summary.log",
    shell:
        "pseudodefense summary --taxon '{params.taxon}' {params.levels} "
        "--assembly-source {params.source} {params.atypical} --out {output} 2> {log}"


rule select_genomes:
    """Reproducible subset (n_genomes, random_seed) of the available assemblies."""
    input:
        rules.genome_summary.output,
    output:
        f"{OUT}/metadata/selected_genomes.tsv",
    params:
        n="" if config["n_genomes"] is None else f"--n {config['n_genomes']}",
        seed=config["random_seed"],
    log:
        f"{OUT}/logs/select_genomes.log",
    shell:
        "pseudodefense select --summary {input} --out {output} {params.n} "
        "--seed {params.seed} > {log} 2>&1"


rule download_genomes:
    """Genome FASTA files plus manifest (accession, strain, level, date, version)."""
    input:
        rules.select_genomes.output,
    output:
        genomes=directory(f"{OUT}/genomes"),
        manifest=f"{OUT}/metadata/manifest.tsv",
    log:
        f"{OUT}/logs/download_genomes.log",
    shell:
        "pseudodefense download --selected {input} --genomes-dir {output.genomes} "
        "--manifest {output.manifest} > {log} 2>&1"
