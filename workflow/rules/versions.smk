rule versions:
    """Record tool versions and the genome download date for this run."""
    input:
        rules.download_genomes.output.manifest,
    output:
        f"{OUT}/versions.txt",
    shell:
        """
        {{
          echo "run_date: $(date -I)"
          echo "genome_download_date: $(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="download_date") c=i}} NR==2{{print $c}}' {input})"
          echo "python: $(python --version 2>&1 | cut -d' ' -f2)"
          echo "snakemake: $(snakemake --version)"
          echo "ncbi-datasets-cli: $(datasets --version | awk '{{print $NF}}')"
          echo "pyrodigal: $(python -c 'import pyrodigal; print(pyrodigal.__version__)')"
          echo "pseudodefense: $(python -c 'import pseudodefense; print(pseudodefense.__version__)')"
        }} > {output}
        """
