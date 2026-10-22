# Cambios

## 0.1.0 (2026-10-22)

Primera version publica.

- Descarga de genomas de NCBI Datasets por lotes, con filtro taxonomico por ANI de NCBI.
- QC de ensamblajes: longitud, secuencias, GC y CheckM.
- Prediccion de genes con pyrodigal y deteccion de sistemas con DefenseFinder 3.0.0, en un entorno conda aislado.
- Parser a tablas tidy, con eliminacion avisada de las filas duplicadas que escribe DefenseFinder.
- Matriz de presencia/ausencia y prevalencia con IC de Wilson.
- Co-ocurrencia con Fisher y Benjamini-Hochberg, con filtro de pares testables y marca de ligamiento fisico (≤ 50 kb).
- Desreplicacion con Mash a 99,9 / 99,8 / 99,5 % ANI y validacion de pares cercanos a los umbrales con FastANI.
- Figuras, informe Markdown/HTML y `versions.txt`.
- Comando `pseudodefense run`, Dockerfile y CI en GitHub Actions (lint, pruebas y dry-run).
