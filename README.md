# pseudodefense

Pipeline reproducible y abierto (Python + Snakemake) para anotar sistemas de defensa antifago en genomas publicos de *Pseudomonas* y describir su prevalencia y co-ocurrencia.

> **Estado:** en desarrollo (hito 4: corrida completa y desreplicacion). No es una novedad cientifica: replica y extiende analisis publicados (ver [docs/REFERENCES.md](docs/REFERENCES.md)).

## Organizacion de carpetas

El repo vive dentro de una carpeta de proyecto que agrupa todo lo que el pipeline genera o necesita y que no va a GitHub:

```
pseudodefense_proyecto/
├── pseudodefense/      este repo (codigo)
│   └── results/        una subcarpeta por corrida: full/, pilot100/, test20/ (no se versiona)
├── miniforge3/         conda (entorno principal `pseudodefense`)
├── conda_envs/         entornos de las reglas (DefenseFinder, Mash), creados por Snakemake
├── resources/          modelos de DefenseFinder
└── logs/               logs de las corridas lanzadas a mano
```

Las rutas a `conda_envs/` y `resources/` son relativas al repo (`../`) y estan en `profiles/default/config.yaml` y `config/config.yaml`.

## Instalacion (Linux / WSL)

```bash
mkdir pseudodefense_proyecto && cd pseudodefense_proyecto
git clone <url-del-repo> pseudodefense
# Miniforge en la carpeta del proyecto: https://github.com/conda-forge/miniforge
bash Miniforge3-Linux-x86_64.sh -b -p "$PWD/miniforge3"
cd pseudodefense
../miniforge3/bin/mamba env create -f environment.yml
source ../miniforge3/bin/activate pseudodefense
```

## Uso

```bash
snakemake -n                 # ver el plan
snakemake --cores 12         # corrida completa en results/full
snakemake --cores 12 --config outdir=results/pilot100 n_genomes=100   # prueba reducida
```

El perfil `profiles/default` se aplica solo: activa los entornos conda de cada regla y los guarda en `../conda_envs`. Toda la configuracion del analisis (taxon, filtros, numero de genomas, semilla, umbrales de QC, desreplicacion) esta en `config/config.yaml`.

## Limitaciones

- DefenseFinder detecta **sistemas conocidos**: la ausencia de deteccion no es ausencia de defensa.
- Los genomas de NCBI sobrerrepresentan aislados clinicos de *P. aeruginosa*; los resultados son descriptivos, no causales.
- La co-ocurrencia (Fisher + Benjamini-Hochberg) trata a los genomas como independientes, aunque no lo son (aislados casi identicos, clados). Los pares de sistemas en el mismo locus se marcan como fisicamente ligados: su co-ocurrencia es esperable.

## Licencia

MIT
