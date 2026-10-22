# pseudodefense: imagen con el entorno principal y los entornos de las reglas ya creados.
# Mantiene la misma organizacion de carpetas que una instalacion local:
#   /opt/pseudodefense_proyecto/{pseudodefense, conda_envs, resources}
# Uso (desde la carpeta donde quieras los resultados):
#   docker build -t pseudodefense .
#   docker run --rm -v "$PWD/results:/opt/pseudodefense_proyecto/pseudodefense/results" \
#     -v "$PWD/resources:/opt/pseudodefense_proyecto/resources" \
#     pseudodefense pseudodefense run --cores 8 --config outdir=results/pilot n_genomes=20
FROM condaforge/miniforge3:24.9.2-0

WORKDIR /opt/pseudodefense_proyecto/pseudodefense
COPY . .

RUN mamba env create -f environment.yml && mamba clean -afy

SHELL ["conda", "run", "--no-capture-output", "-n", "pseudodefense", "/bin/bash", "-c"]

# Entornos de DefenseFinder, Mash y FastANI en ../conda_envs (profiles/default).
# Los modelos de DefenseFinder se descargan en la primera corrida a ../resources
# (montar ese volumen los conserva entre corridas; su version queda en versions.txt).
RUN snakemake --conda-create-envs-only --cores 1 && mamba clean -afy

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "pseudodefense"]
CMD ["pseudodefense", "run", "--cores", "4"]
