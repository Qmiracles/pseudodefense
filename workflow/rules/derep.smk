DEREP = config["derep"]
DEREP_LABELS = [f"ani{t:g}" for t in DEREP["thresholds"]]


rule mash_sketch:
    """Mash sketches of every genome that passed QC."""
    input:
        accessions=f"{OUT}/qc/passed_accessions.txt",
        genomes=rules.download_genomes.output.genomes,
    output:
        sketch=f"{OUT}/derep/genomes.msh",
        paths=f"{OUT}/derep/genome_paths.txt",
    params:
        s=DEREP["sketch_size"],
        k=DEREP["kmer"],
    threads: config["threads"]
    conda:
        "../envs/mash.yml"
    log:
        f"{OUT}/logs/mash_sketch.log",
    shell:
        "sed 's|^|{input.genomes}/|; s|$|.fna.gz|' {input.accessions} > {output.paths} && "
        "mash sketch -s {params.s} -k {params.k} -p {threads} -l {output.paths} "
        "-o {OUT}/derep/genomes > {log} 2>&1"


rule mash_dist:
    """All-vs-all Mash distances (gzipped)."""
    input:
        rules.mash_sketch.output.sketch,
    output:
        f"{OUT}/derep/mash_dist.tsv.gz",
    threads: config["threads"]
    conda:
        "../envs/mash.yml"
    log:
        f"{OUT}/logs/mash_dist.log",
    shell:
        "mash dist -p {threads} {input} {input} 2> {log} | gzip > {output}"


rule dereplicate:
    """Clusters at each ANI threshold and one representative per cluster."""
    input:
        dist=rules.mash_dist.output,
        manifest=rules.download_genomes.output.manifest,
        qc=f"{OUT}/qc/qc_report.tsv",
    output:
        clusters=f"{OUT}/derep/clusters.tsv",
        summary=f"{OUT}/derep/summary.tsv",
        reps=expand(f"{OUT}/derep/representatives_{{label}}.txt", label=DEREP_LABELS),
    params:
        thresholds=" ".join(f"--threshold {t}" for t in DEREP["thresholds"]),
    log:
        f"{OUT}/logs/dereplicate.log",
    shell:
        "pseudodefense dereplicate --dist {input.dist} --manifest {input.manifest} "
        "--qc-report {input.qc} {params.thresholds} --outdir {OUT}/derep > {log} 2>&1"
