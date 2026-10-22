DF_MODELS = config["defensefinder"]["models_dir"]


def defensefinder_outputs(wildcards):
    """DefenseFinder systems tables for every genome that passed QC."""
    return expand(
        f"{OUT}/defensefinder/{{accession}}/{{accession}}_defense_finder_systems.tsv",
        accession=passed_accessions(wildcards),
    )


def gff_outputs(wildcards):
    """pyrodigal GFF3 for every genome that passed QC (needed for coordinates)."""
    return expand(f"{OUT}/genes/{{accession}}.gff", accession=passed_accessions(wildcards))


rule parse_defense:
    """Tidy tables: one row per system and one row per gene, all genomes."""
    input:
        accessions=f"{OUT}/qc/passed_accessions.txt",
        systems=defensefinder_outputs,
        gffs=gff_outputs,
    output:
        systems=f"{OUT}/tables/defense_systems.tsv",
        genes=f"{OUT}/tables/defense_genes.tsv",
    log:
        f"{OUT}/logs/parse_defense.log",
    shell:
        "pseudodefense parse-defense --accessions {input.accessions} "
        "--defensefinder-dir {OUT}/defensefinder --genes-dir {OUT}/genes "
        "--systems-out {output.systems} --genes-out {output.genes} > {log} 2>&1"


rule defensefinder_models:
    """Download DefenseFinder models into the project (versioned per run)."""
    output:
        directory(DF_MODELS),
    conda:
        "../envs/defensefinder.yml"
    log:
        f"{OUT}/logs/defensefinder_models.log",
    shell:
        "defense-finder update --models-dir {output} > {log} 2>&1"


rule defensefinder_versions:
    """Versions of DefenseFinder, its engine and the models used in this run."""
    input:
        DF_MODELS,
    output:
        f"{OUT}/defensefinder_versions.txt",
    conda:
        "../envs/defensefinder.yml"
    shell:
        """
        {{
          echo "defense-finder: $(python -c 'from importlib.metadata import version; print(version("mdmparis-defense-finder"))')"
          echo "macsyfinder: $(python -c 'from importlib.metadata import version; print(version("macsyfinder"))')"
          echo "hmmer: $(hmmsearch -h | grep -m1 -o 'HMMER [0-9.]*' | cut -d' ' -f2)"
          for m in {input}/*/metadata.yml; do
            echo "models_$(basename $(dirname $m)): $(awk '/^vers:/{{print $2}}' $m)"
          done
          echo "models_download_date: $(date -I -r {input})"
        }} > {output}
        """


rule defensefinder:
    """Run DefenseFinder on one proteome."""
    input:
        faa=rules.predict_genes.output.faa,
        models=DF_MODELS,
    output:
        systems=f"{OUT}/defensefinder/{{accession}}/{{accession}}_defense_finder_systems.tsv",
        genes=f"{OUT}/defensefinder/{{accession}}/{{accession}}_defense_finder_genes.tsv",
        hmmer=f"{OUT}/defensefinder/{{accession}}/{{accession}}_defense_finder_hmmer.tsv",
    params:
        outdir=f"{OUT}/defensefinder/{{accession}}",
        db_type=config["defensefinder"]["db_type"],
        coverage=config["defensefinder"]["coverage"],
    threads: config["defensefinder"]["workers"]
    conda:
        "../envs/defensefinder.yml"
    log:
        f"{OUT}/logs/defensefinder/{{accession}}.log",
    benchmark:
        f"{OUT}/benchmarks/defensefinder/{{accession}}.tsv"
    shell:
        "defense-finder run {input.faa} --out-dir {params.outdir} --workers {threads} "
        "--db-type {params.db_type} --coverage {params.coverage} "
        "--models-dir {input.models} > {log} 2>&1"
