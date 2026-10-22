from snakemake.utils import min_version

min_version("9.0")


configfile: "config/config.yaml"


OUT = config["outdir"]


include: "workflow/rules/download.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/genes.smk"
include: "workflow/rules/defense.smk"
include: "workflow/rules/versions.smk"


rule all:
    input:
        f"{OUT}/metadata/manifest.tsv",
        f"{OUT}/qc/qc_report.tsv",
        f"{OUT}/tables/defense_systems.tsv",
        f"{OUT}/versions.txt",
