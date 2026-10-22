# pseudodefense

Pipeline reproducible y abierto (Python + Snakemake) para anotar sistemas de defensa antifago en genomas publicos de *Pseudomonas* y describir su prevalencia y co-ocurrencia.

> **Estado:** en desarrollo (hito 1: descarga + QC). No es una novedad cientifica: replica y extiende analisis publicados (ver [docs/REFERENCES.md](docs/REFERENCES.md)).

## Instalacion (Linux / WSL)

```bash
mamba env create -f environment.yml
conda activate pseudodefense
```

## Uso

```bash
snakemake -n                 # ver el plan
snakemake --cores 4          # descarga + QC de los genomas definidos en config/config.yaml
```

Toda la configuracion (taxon, filtros, numero de genomas, semilla, umbrales de QC) esta en `config/config.yaml`.

## Limitaciones

- DefenseFinder detecta **sistemas conocidos**: la ausencia de deteccion no es ausencia de defensa.
- Los genomas de NCBI sobrerrepresentan aislados clinicos de *P. aeruginosa*; los resultados son descriptivos, no causales.
- La co-ocurrencia (Fisher + Benjamini-Hochberg) trata a los genomas como independientes, aunque no lo son (aislados casi identicos, clados). Los pares de sistemas en el mismo locus se marcan como fisicamente ligados: su co-ocurrencia es esperable.

## Licencia

MIT
