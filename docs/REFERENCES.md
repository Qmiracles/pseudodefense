# Referencias y decisiones metodologicas

Formato: **Decision | Justificacion | Referencia (DOI/URL)**.
Regla: ninguna referencia se agrega de memoria. Lo no verificado se marca `[POR VERIFICAR]`.

## Herramientas

| Herramienta | Uso | Referencia |
|---|---|---|
| DefenseFinder 3.0.0 (modelos defense-finder-models 3.1.0 y CasFinder 3.1.0) | Deteccion de sistemas de defensa (primaria) | Tesson F. et al. (2022). Nat Commun 13. https://doi.org/10.1038/s41467-022-30269-9 |
| PADLOC (CLI) | Comparacion entre herramientas (opcional) | Payne L.J. et al. (2021). Nucleic Acids Res 49:10868-10878. https://doi.org/10.1093/nar/gkab883 |
| PADLOC (web / RefSeq precalculado) | Solo si se usan sus resultados precalculados | Payne L.J. et al. (2022). Nucleic Acids Res 50:W541-W550. https://doi.org/10.1093/nar/gkac400 |
| Prodigal | Algoritmo de prediccion de genes | Hyatt D. et al. (2010). BMC Bioinformatics 11:119. https://doi.org/10.1186/1471-2105-11-119 |
| pyrodigal 3.7 | Binding de Python a Prodigal (citar ambos) | Larralde M. (2022). JOSS 7(72):4296. https://doi.org/10.21105/joss.04296 |
| FastANI | Desreplicacion por ANI (opcional) | Jain C. et al. (2018). Nat Commun 9:5114. https://doi.org/10.1038/s41467-018-07641-9 |
| Snakemake 9.27 | Orquestacion | Molder F. et al. (2021). F1000Research 10:33. https://f1000research.com/articles/10-33 |
| BACPHLIP (opcional) | Estilo de vida de fagos | Hockenberry A.J. & Wilke C.O. (2021). PeerJ 9:e11396. https://doi.org/10.7717/peerj.11396 |
| NCBI Datasets CLI 18.39 | Descarga de genomas | Documentacion: https://github.com/ncbi/datasets ; version registrada en `manifest.tsv` y `versions.txt`. `[POR VERIFICAR si existe articulo citable]` |
| MacSyFinder | Motor de DefenseFinder | `[POR VERIFICAR autores y ano]` |

## Decisiones

| Decision | Justificacion | Referencia |
|---|---|---|
| DefenseFinder como herramienta primaria | Detecta sistemas conocidos con modelos actualizables; aplicado a >20 000 genomas | Tesson et al. 2022 |
| Reportar que la deteccion es conservadora | La configuracion prioriza sistemas completos; posible subdeteccion | Tesson et al. 2022 |
| DefenseFinder en un entorno conda aislado (`workflow/envs/defensefinder.yml`) | v3.0.0 fija `click==8.0.3` y `macsyfinder==2.1.4` (requires_dist en PyPI, consultado 2026-10-21); aislarlo evita conflictos con typer/snakemake | https://pypi.org/project/mdmparis-defense-finder/3.0.0/ |
| Prediccion de genes con pyrodigal | Resultados identicos a Prodigal v2.6.3; instalable sin binario externo | Larralde 2022; Hyatt 2010 |
| pyrodigal en modo single, entrenado una vez por genoma con todas sus secuencias (cromosoma + plasmidos), tabla 11, sin `*` de stop | Asi lo hace Prodigal con varios contigs (los une con enlazadores para entrenar, segun la documentacion de `GeneFinder.train`); genomas completos de >6 Mb superan el minimo recomendado. El modulo exporta FASTA de proteinas y GFF3 (requisito de PADLOC) | Docstring de pyrodigal 3.7.1 `GeneFinder.train`; Payne et al. 2021 |
| DefenseFinder con `--db-type gembase` sobre las proteinas de pyrodigal | Cada contig se trata como replicon distinto (evita sistemas "pegados" entre el final del cromosoma y un plasmido). DefenseFinder 3.0.0 hace exactamente esto cuando recibe nucleotidos: predice con pyrodigal, escribe IDs `<contig>_<n>` y fija `dbtype = "gembase"` | Codigo fuente `defense_finder_cli/main.py` (v3.0.0), lineas 186-219 |
| Modelos de DefenseFinder descargados en la carpeta del proyecto (`../resources/defensefinder_models`, fuera del repo) y version registrada en `versions.txt` | Reproducibilidad y comparabilidad con otros trabajos; Georjon et al. reportan version de herramienta y modelos | Georjon et al. 2023. https://doi.org/10.1099/mic.0.001374 |
| Parametros de DefenseFinder por defecto (cobertura 0.4, umbrales GA de los perfiles) | No se modifican los umbrales validados por los autores | Tesson et al. 2022; `defense-finder run --help` |
| Eliminar filas exactamente duplicadas de `*_defense_finder_systems.tsv` (con aviso en el log) | DefenseFinder 3.0.0 escribio la misma linea dos veces en al menos un genoma (GCF_012971705.1, sistema `DNA_met_RM_Type_II_PDDEXK_13`); filas identicas no aportan informacion. Filas distintas con el mismo `sys_id` detienen el parser | Observacion directa en este proyecto (fixture en `tests/data/`) |
| Conservar proteinas compartidas entre sistemas | DefenseFinder puede asignar una proteina a dos sistemas solapados (p. ej. una REasa en dos RM tipo II); se respeta su asignacion y queda visible en la tabla de genes | Observacion directa; columna `used_in` de DefenseFinder |
| Genomas RefSeq, nivel `complete` + `chromosome`, `--exclude-atypical` | Ensamblajes cerrados reducen artefactos en la deteccion de sistemas multigenicos; RefSeq aporta metadatos curados (incluido CheckM) | Decision del proyecto. Flags: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/ |
| Muestreo reproducible de `n_genomes` (semilla fija, orden por accesion) | Independiza la seleccion del orden de respuesta de NCBI | Decision del proyecto |
| QC: longitud 5.5-7.8 Mb, <=20 secuencias, GC 65-67.5 % | Umbrales provisionales para *P. aeruginosa*; revisar con la distribucion observada del dataset | `[POR VERIFICAR: buscar referencia de tamano y GC de P. aeruginosa]` |
| QC: CheckM completitud >=95 %, contaminacion <=5 % (si NCBI lo reporta) | Umbrales de uso comun para genomas de alta calidad | `[POR VERIFICAR: estandar MIMAG / referencia citable]` |
| Presencia/ausencia a nivel de `type` de DefenseFinder (configurable a `subtype`); conteos por genoma incluyen genomas sin sistemas | El nivel `type` agrupa variantes del mismo mecanismo; mantener genomas con 0 sistemas evita sesgar prevalencias al alza | Decision del proyecto |
| Intervalo de Wilson al 95 % para prevalencias | Mejor cobertura que el intervalo de Wald con n pequeno o proporciones cercanas a 0 o 1 | `[POR VERIFICAR: Wilson 1927 J Am Stat Assoc 22:209-212; Brown, Cai & DasGupta 2001 Stat Sci 16:101-133]` |
| Co-ocurrencia: prueba exacta de Fisher bilateral por par de sistemas | Exacta con conteos pequenos; no asume aproximacion asintotica | `[POR VERIFICAR referencia citable]`; implementacion `scipy.stats.fisher_exact` |
| Solo se prueban sistemas con prevalencia 5-95 % y presentes y ausentes en >= 3 genomas | Pares con sistemas casi universales o muy raros no pueden alcanzar significancia pero inflan el numero de pruebas en BH. Con los 20 genomas de prueba: 71 sistemas / 2 485 pares sin el minimo absoluto; 26 / 325 con el | Decision del proyecto; idea relacionada con pruebas "no testables" `[POR VERIFICAR: Tarone 1990, Biometrics]` |
| log2 OR con correccion de Haldane-Anscombe (+0.5) solo para graficar | El OR muestral es 0 o infinito con celdas vacias. Limitacion: con conteos pequenos el estimador corregido esta sesgado y puede tener signo opuesto al OR muestral; la direccion de la asociacion se toma del OR muestral | `[POR VERIFICAR referencia citable]` |
| Marcar pares fisicamente ligados (distancia mediana <= 10 kb en el mismo replicon) | Sistemas en el mismo locus co-ocurren trivialmente. En los 20 genomas de prueba, los 2 pares significativos (JukAB-Ophion, separados ~250 pb; CoCoNut-GmrSD_RM_Type_IV, solapados) son ligados. Umbral de 10 kb: decision del proyecto, revisar con la corrida completa | Observacion directa; umbral `[POR DEFINIR]` |
| No se corrige por estructura poblacional/filogenia en las pruebas de co-ocurrencia | Los genomas no son independientes (aislados casi identicos, clados); los p-valores sobreestiman la evidencia. Mitigacion prevista: desreplicacion (hito 4) y reportar resultados como descriptivos | Limitacion declarada |
| Umbral de 95 % ANI para especie | Discontinuidad intra (>95 %) vs inter-especie (<83 %) | Jain et al. 2018 |
| Umbral de desreplicacion | Sin cifra citable universal; analisis de sensibilidad con varios umbrales | `[POR DEFINIR]` |
| Benjamini-Hochberg en co-ocurrencia | Controla FDR; asume independencia (limitacion a mencionar) | Benjamini & Hochberg 1995. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x |
| Usar Snakemake | Reproducibilidad y transparencia | Molder et al. 2021 |
| (Opcional) Solo genomas completos para BACPHLIP | Con genomas incompletos clasifica como virulento por defecto | Hockenberry & Wilke 2021; https://pypi.org/project/bacphlip |

## Literatura de contexto sobre *Pseudomonas*

1. Core defense hotspots within *P. aeruginosa*... Nucleic Acids Res 2023. https://doi.org/10.1093/nar/gkad317
2. Accumulation of defense systems in phage-resistant strains of *P. aeruginosa*. https://pmc.ncbi.nlm.nih.gov/articles/PMC10889362 `[POR VERIFICAR autores y revista]`
3. Defence systems drive accessory genome interactions in *P. aeruginosa*. ISME Communications. https://doi.org/10.1093/ismeco/ycag130
4. Co-occurrence analysis of anti-phage systems in *P. aeruginosa* reveals prevalent neutral interactions. Microorganisms 13, 2783 (2025). `[POR VERIFICAR autores]`
5. Bacteriophage genomes encode both broad and specific counter-defense repertoires... https://repository.tudelft.nl/file/File_7e1db523-5d5f-4330-9d17-f1c4315a03fb `[POR VERIFICAR autores y revista]`
6. Georjon H. et al. (2023). Genomic characterization of the antiviral arsenal of Actinobacteria. Microbiology. https://doi.org/10.1099/mic.0.001374
