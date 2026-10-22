rule predict_genes:
    """Proteins (.faa, for DefenseFinder) and GFF3 (for PADLOC) with pyrodigal."""
    input:
        f"{OUT}/genomes/{{accession}}.fna",
    output:
        faa=f"{OUT}/genes/{{accession}}.faa",
        gff=f"{OUT}/genes/{{accession}}.gff",
    params:
        table=config["genes"]["translation_table"],
    log:
        f"{OUT}/logs/genes/{{accession}}.log",
    shell:
        "pseudodefense genes --fasta {input} --faa {output.faa} --gff {output.gff} "
        "--translation-table {params.table} > {log} 2>&1"
