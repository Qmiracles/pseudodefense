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


rule mash_version:
    """Mash version from the rule environment (for versions.txt)."""
    output:
        f"{OUT}/derep/mash_version.txt",
    conda:
        "../envs/mash.yml"
    shell:
        'echo "mash: $(mash --version)" > {output}'


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


VALID = DEREP["validation"]


rule ani_validation_pairs:
    """Sample genome pairs whose Mash ANI is near each dereplication threshold."""
    input:
        dist=rules.mash_dist.output,
    output:
        f"{OUT}/derep/ani_validation_pairs.tsv",
    params:
        thresholds=" ".join(f"--threshold {t}" for t in DEREP["thresholds"]),
        n=VALID["pairs_per_threshold"],
        window=VALID["window"],
        seed=config["random_seed"],
    log:
        f"{OUT}/logs/ani_validation_pairs.log",
    shell:
        "pseudodefense ani-pairs --dist {input.dist} {params.thresholds} "
        "--per-threshold {params.n} --window {params.window} --seed {params.seed} "
        "--out {output} > {log} 2>&1"


rule fastani:
    """Alignment-based ANI (FastANI) for the validation pairs."""
    input:
        pairs=rules.ani_validation_pairs.output,
        genomes=rules.download_genomes.output.genomes,
    output:
        raw=f"{OUT}/derep/fastani_raw.tsv",
        version=f"{OUT}/derep/fastani_version.txt",
    threads: config["threads"]
    conda:
        "../envs/fastani.yml"
    log:
        f"{OUT}/logs/fastani.log",
    shell:
        """
        echo "fastani: $(fastANI --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1)" > {output.version}
        tmp=$(mktemp -d {OUT}/derep/fastani_tmp.XXXX)
        tail -n +2 {input.pairs} | while IFS=$'\\t' read -r a b rest; do
          rm -f $tmp/out.tsv
          fastANI -q {input.genomes}/$a.fna.gz -r {input.genomes}/$b.fna.gz \
            -t {threads} -o $tmp/out.tsv >> {log} 2>&1
          if [ -s $tmp/out.tsv ]; then cat $tmp/out.tsv; fi
        done > {output.raw}
        rm -r $tmp
        """


rule ani_validation:
    """Compare Mash and FastANI near the thresholds."""
    input:
        pairs=rules.ani_validation_pairs.output,
        fastani=rules.fastani.output.raw,
    output:
        table=f"{OUT}/derep/ani_validation.tsv",
        summary=f"{OUT}/derep/ani_validation_summary.tsv",
    log:
        f"{OUT}/logs/ani_validation.log",
    shell:
        "pseudodefense ani-compare --pairs {input.pairs} --fastani {input.fastani} "
        "--table-out {output.table} --summary-out {output.summary} > {log} 2>&1"
