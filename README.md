# pseudodefense

[![CI](https://github.com/Qmiracles/pseudodefense/actions/workflows/ci.yml/badge.svg)](https://github.com/Qmiracles/pseudodefense/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23577682.svg)](https://doi.org/10.5281/zenodo.23577682)

Pipeline reproducible y abierto (Python + Snakemake) para anotar sistemas de defensa antifago en genomas publicos de *Pseudomonas* y describir su prevalencia y co-ocurrencia.

> **Estado:** v0.1.0. No es una novedad cientifica: aplica herramientas publicadas (DefenseFinder, Mash, FastANI) a datos publicos de NCBI, y su aporte es que el analisis es reproducible y actualizable. Los trabajos cercanos estan en [docs/REFERENCES.md](docs/REFERENCES.md).

## Que hace

Un solo comando ejecuta, de principio a fin:

1. **Descarga** de NCBI los ensamblajes de un taxon (por defecto *P. aeruginosa*, RefSeq, nivel `complete` y `chromosome`) que pasan el control taxonomico por ANI de NCBI.
2. **Control de calidad**: longitud, numero de secuencias, GC y CheckM.
3. **Prediccion de genes** con pyrodigal.
4. **Deteccion de sistemas de defensa** con DefenseFinder, un genoma por trabajo y en paralelo.
5. **Tablas tidy** de sistemas y genes, y una matriz de presencia/ausencia.
6. **Desreplicacion** con Mash a varios umbrales de ANI (99,9 / 99,8 / 99,5 %). Una muestra de pares cercanos a cada umbral se re-mide con FastANI: en la corrida de referencia FastANI dio valores unos 0,07 puntos mas altos, asi que los umbrales estan en la escala de Mash.
7. **Estadistica**:
   - prevalencia con intervalo de Wilson;
   - co-ocurrencia con Fisher y Benjamini-Hochberg, con todos los genomas y con cada conjunto desreplicado;
   - marca de los pares de sistemas fisicamente ligados (≤ 50 kb).
8. **Figuras** (PNG + SVG), un **informe** (`report.md` y `report.html`) y `versions.txt` con las versiones de todas las herramientas y modelos.

## Resultados de la corrida de referencia

Corrida del 22 de octubre de 2026 sobre 1 847 genomas de *P. aeruginosa* (DefenseFinder 3.0.0, modelos 3.1.0):

- **Repertorio:** 25 093 sistemas de 181 tipos; mediana de 12 sistemas por genoma (rango 2-33).
- **Tipos mas frecuentes:**
  - RM tipo I: 59 %;
  - CRISPR-Cas: 48 % (sobre todo I-F);
  - Gabija: 40 %;
  - CBASS: 38 %.
- **Casi universales:** Lugos y Dag aparecen en el 99 % de los genomas. Estan pendientes de verificar: pueden ser genes del genoma central.
- **Efecto de la desreplicacion:**
  - con un genoma por cluster de 99,5 % ANI quedan 362 genomas;
  - los pares de sistemas asociados bajan de 663 a 54;
  - las exclusiones desaparecen.
- **Ligamiento fisico:** 36 de esos 54 pares tienen los dos sistemas a ≤ 50 kb, sobre todo un bloque Ophion-Wadjet-JukAB-CoCoNut-GmrSD-Druantia.

Los resultados (12 GB) no se versionan; se regeneran con el pipeline.

## Instalacion (Linux / WSL)

Todo vive en una carpeta de proyecto; el repo es una subcarpeta.

```bash
mkdir pseudodefense_proyecto && cd pseudodefense_proyecto
git clone https://github.com/Qmiracles/pseudodefense.git
# Miniforge dentro de la carpeta del proyecto: https://github.com/conda-forge/miniforge
bash Miniforge3-Linux-x86_64.sh -b -p "$PWD/miniforge3"
cd pseudodefense
../miniforge3/bin/mamba env create -f environment.yml
source ../miniforge3/bin/activate pseudodefense
```

Los entornos de DefenseFinder, Mash y FastANI los crea Snakemake en la primera corrida, en `../conda_envs`. DefenseFinder va aparte porque fija versiones de dependencias incompatibles con el entorno principal. Los modelos de DefenseFinder se descargan en `../resources`.

```
pseudodefense_proyecto/
├── pseudodefense/      este repo
│   └── results/        una subcarpeta por corrida (no se versiona)
├── miniforge3/         conda (entorno principal `pseudodefense`)
├── conda_envs/         entornos de las reglas, creados por Snakemake
└── resources/          modelos de DefenseFinder
```

### Docker

```bash
docker build -t pseudodefense .
docker run --rm \
  -v "$PWD/results:/opt/pseudodefense_proyecto/pseudodefense/results" \
  -v "$PWD/resources:/opt/pseudodefense_proyecto/resources" \
  pseudodefense pseudodefense run --cores 8 --config outdir=results/pilot n_genomes=20
```

## Uso

```bash
pseudodefense run -n                    # ver el plan sin ejecutar nada
pseudodefense run --cores 12            # corrida completa en results/full
pseudodefense run --cores 12 --config outdir=results/pilot100 n_genomes=100   # prueba reducida
```

`pseudodefense run` llama a Snakemake con `--rerun-incomplete`. Si la corrida se interrumpe, al relanzarla retoma desde donde quedo. Equivale a `snakemake --cores N` desde la raiz del repo.

Requisitos medidos en la corrida completa:

- 1 847 genomas en 10 h 34 min con 12 nucleos y 15 GB de RAM;
- 12 GB de disco;
- unos 85 s y 180 MB de RAM por genoma en DefenseFinder.

`scripts/run_full.sh` la lanza desacoplada de la terminal.

Toda la configuracion vive en [`config/config.yaml`](config/config.yaml):

- taxon y filtros de NCBI;
- numero de genomas y semilla;
- umbrales de QC;
- nivel de analisis (tipo o subtipo);
- umbrales de co-ocurrencia y de ligamiento;
- umbrales de desreplicacion y de su validacion;
- genomas de referencia del informe.

### Salidas principales (`results/<corrida>/`)

| Archivo | Contenido |
| --- | --- |
| `report.html`, `report.md` | Informe de la corrida |
| `metadata/manifest.tsv` | Un genoma por fila: acceso, cepa, nivel, CheckM, ANI de NCBI, fecha de descarga |
| `qc/qc_report.tsv` | Metricas y motivo de descarte de cada genoma |
| `tables/defense_systems.tsv`, `tables/defense_genes.tsv` | Sistemas y genes detectados, con coordenadas |
| `tables/matrix_type.tsv` | Matriz genomas x tipos (0/1) |
| `stats/prevalence_by_genome_set.tsv` | Prevalencia (IC de Wilson) con todos los genomas y por umbral de ANI |
| `stats/cooccurrence*.tsv` | Fisher, q de BH, distancia entre sistemas y ligamiento |
| `derep/clusters.tsv`, `derep/ani_validation*.tsv` | Clusters, representantes y validacion Mash vs FastANI |
| `figures/` | Mapa de presencia, prevalencia, sistemas por genoma, co-ocurrencia |
| `versions.txt` | Fechas y versiones de herramientas y modelos |

## Desarrollo

```bash
pytest -q                        # pruebas
ruff check src tests && ruff format --check src tests
snakemake -n --workflow-profile none --config outdir=results/ci n_genomes=2   # lo mismo que el CI
```

## Limitaciones

- **Solo sistemas conocidos.** DefenseFinder detecta sistemas con modelo publicado: la ausencia de deteccion no es ausencia de defensa.
- **Sesgo de muestreo.** NCBI sobrerrepresenta aislados clinicos y linajes muy secuenciados. La desreplicacion lo reduce, pero los resultados son descriptivos, no causales.
- **Genomas no independientes.** La co-ocurrencia no corrige por filogenia; los pares fisicamente ligados se marcan porque su co-ocurrencia es esperable.
- **Umbrales provisionales.** Los de QC estan pendientes de referencia citable. Los de desreplicacion y ligamiento salen de los propios datos (ver [docs/REFERENCES.md](docs/REFERENCES.md)).

## Como citar

Quevedo, A. (2026). *pseudodefense: pipeline reproducible para sistemas de defensa antifago en genomas de Pseudomonas* (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.23577683

El DOI de concepto [10.5281/zenodo.23577682](https://doi.org/10.5281/zenodo.23577682) apunta siempre a la ultima version. Metadatos en [`CITATION.cff`](CITATION.cff). Cita tambien las herramientas que usa el pipeline (DefenseFinder, Prodigal/pyrodigal, Mash, FastANI, Snakemake), listadas en [docs/REFERENCES.md](docs/REFERENCES.md).

## Licencia

[MIT](LICENSE)
