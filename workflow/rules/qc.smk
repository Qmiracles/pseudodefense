checkpoint qc:
    """Assembly QC. Checkpoint: later per-genome rules use the accessions that pass."""
    input:
        manifest=rules.download_genomes.output.manifest,
        genomes=rules.download_genomes.output.genomes,
    output:
        report=f"{OUT}/qc/qc_report.tsv",
        passed=f"{OUT}/qc/passed_accessions.txt",
        summary=f"{OUT}/qc/qc_summary.txt",
    params:
        args=" ".join(f"--{k.replace('_', '-')} {v}" for k, v in config["qc"].items()),
    log:
        f"{OUT}/logs/qc.log",
    shell:
        "pseudodefense qc --manifest {input.manifest} --report {output.report} "
        "--passed {output.passed} {params.args} > {output.summary} 2> {log}"


def passed_accessions(wildcards):
    """Accessions that passed QC (resolved after the checkpoint runs)."""
    path = checkpoints.qc.get().output.passed
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]
