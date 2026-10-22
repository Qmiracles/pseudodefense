"""Run report: one Markdown and one self-contained HTML file from the pipeline outputs.

The report is a list of blocks (heading, paragraph, table, image, preformatted
text) rendered twice, so both formats always say the same thing. It only reads
files the pipeline already wrote; it computes no new statistics.
"""

from __future__ import annotations

import base64
import html
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

NEAR_UNIVERSAL = 0.95  # prevalence above which a type is flagged for review


@dataclass
class Block:
    """One element of the report."""

    kind: str  # "h1" | "h2" | "p" | "table" | "img" | "pre"
    text: str = ""
    table: pd.DataFrame | None = None


def _fmt(value: object) -> str:
    """Short text for a table cell: integers plain, numbers below 0.1 in 2 significant digits.

    Two significant digits keep a q-value of 0.0477 from printing as 0.05.
    """
    if isinstance(value, float):
        if pd.isna(value):
            return ""
        if value.is_integer():
            return str(int(value))
        if abs(value) < 0.1:
            return f"{value:.2g}"
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return str(value)


def md_table(df: pd.DataFrame) -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table."""
    header = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join(_fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([header, rule, *rows])


def _read_versions(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            out[key.strip()] = value.strip()
    return out


def _pct(x: float) -> float:
    return round(100 * x, 1)


def _set_label(name: str) -> str:
    return "todos" if name == "all" else name.replace("ani", "") + " % ANI (Mash)"


def build_report(
    results: Path,
    taxon: str,
    cooc_set: str,
    references: dict[str, str],
    alpha: float = 0.05,
    linkage_bp: int = 50_000,
) -> list[Block]:
    """Assemble the report blocks from a run directory (e.g. ``results/full``)."""
    versions = _read_versions(results / "versions.txt")
    manifest = pd.read_csv(results / "metadata/manifest.tsv", sep="\t")
    qc = pd.read_csv(results / "qc/qc_report.tsv", sep="\t")
    systems = pd.read_csv(results / "tables/defense_systems.tsv", sep="\t")
    per_genome = pd.read_csv(results / "tables/systems_per_genome.tsv", sep="\t")
    prev = pd.read_csv(results / "stats/prevalence_by_genome_set.tsv", sep="\t")
    derep_summary = pd.read_csv(results / "derep/summary.tsv", sep="\t")
    clusters = pd.read_csv(results / "derep/clusters.tsv", sep="\t")
    ani_summary_path = results / "derep/ani_validation_summary.tsv"

    blocks = [
        Block("h1", f"pseudodefense: informe de la corrida ({taxon})"),
        Block(
            "p",
            f"Corrida del {versions.get('run_date', '?')}; genomas descargados de NCBI el "
            f"{versions.get('genome_download_date', '?')}; pseudodefense "
            f"{versions.get('pseudodefense', '?')}, DefenseFinder "
            f"{versions.get('defense-finder', '?')} con modelos "
            f"{versions.get('models_defense-finder-models', '?')}. Este informe solo resume "
            "archivos que el pipeline ya escribio en esta carpeta.",
        ),
    ]

    # --- Data and QC
    passed = qc["passed"].sum()
    levels = manifest["assembly_level"].value_counts()
    blocks += [
        Block("h2", "Datos y control de calidad"),
        Block(
            "p",
            f"{len(qc)} genomas evaluados, {passed} aprobados y {len(qc) - passed} descartados. "
            "Niveles de ensamblaje descargados: "
            + ", ".join(f"{n} {lvl}" for lvl, n in levels.items())
            + ".",
        ),
    ]
    failed = qc.loc[~qc["passed"], ["accession", "reasons"]]
    if len(failed):
        blocks.append(Block("table", table=failed.rename(columns={
            "accession": "Genoma", "reasons": "Motivo"})))  # fmt: skip

    # --- Systems
    n = per_genome["n_systems"]
    prev_all = prev[prev["genome_set"] == "all"].sort_values("prevalence", ascending=False)
    universal = prev_all[prev_all["prevalence"] > NEAR_UNIVERSAL]
    blocks += [
        Block("h2", "Sistemas de defensa"),
        Block(
            "p",
            f"{len(systems)} sistemas en {systems['accession'].nunique()} genomas: "
            f"{systems['type'].nunique()} tipos y {systems['subtype'].nunique()} subtipos. "
            f"Sistemas por genoma: mediana {n.median():g} (cuartiles {n.quantile(0.25):g} y "
            f"{n.quantile(0.75):g}; rango {n.min()}-{n.max()}).",
        ),
    ]
    if len(universal):
        names = ", ".join(
            f"{r.system} ({_pct(r.prevalence)} %)" for r in universal.itertuples(index=False)
        )
        blocks.append(
            Block(
                "p",
                f"Tipos presentes en mas del {NEAR_UNIVERSAL:.0%} de los genomas: {names}. "
                "Un sistema casi universal puede ser una defensa conservada o un modelo que "
                "reconoce un gen del genoma central; conviene revisarlo antes de interpretarlo. "
                "Se excluyen de las pruebas de co-ocurrencia.",
            )
        )

    # --- Prevalence by genome set
    wide = prev.pivot(index="system", columns="genome_set", values="prevalence") * 100
    sets = ["all", *sorted(c for c in wide.columns if c != "all")[::-1]]
    totals = prev.groupby("genome_set")["total_genomes"].first()
    top = wide.loc[prev_all["system"].head(20), sets].round(1).reset_index()
    top.columns = ["Tipo", *[f"{_set_label(s)} (n = {totals[s]})" for s in sets]]
    blocks += [
        Block("h2", "Prevalencia (% de genomas)"),
        Block(
            "p",
            "Los 20 tipos mas frecuentes con todos los genomas y con un genoma por cluster a "
            "cada umbral de ANI. Intervalos de Wilson en stats/prevalence_by_genome_set.tsv.",
        ),
        Block("table", table=top),
    ]

    # --- Dereplication
    ds = derep_summary.rename(columns={
        "threshold_ani": "Umbral ANI Mash (%)", "n_genomes": "Genomas", "n_clusters": "Clusters",
        "largest_cluster": "Cluster mas grande", "singletons": "Clusters de 1 genoma"})  # fmt: skip
    blocks += [Block("h2", "Desreplicacion (Mash)"), Block("table", table=ds)]
    labels = [c for c in clusters.columns if c.startswith("ani")]
    ref_rows = []
    for name, acc in references.items():
        row = clusters[clusters["accession"] == acc]
        if row.empty:
            ref_rows.append({"Genoma": name, "Acceso": acc, "En la corrida": "no"})
            continue
        r = row.iloc[0]
        entry = {"Genoma": name, "Acceso": acc, "En la corrida": "si"}
        for lab in labels:
            size = int((clusters[lab] == r[lab]).sum())
            rep = "representante" if r[f"rep_{lab}"] else "no representante"
            entry[_set_label(lab)] = f"cluster de {size} ({rep})"
        ref_rows.append(entry)
    if ref_rows:
        blocks += [
            Block("p", "Ubicacion de los genomas de referencia en los clusters:"),
            Block("table", table=pd.DataFrame(ref_rows).fillna("")),
        ]

    # --- ANI validation
    if ani_summary_path.exists():
        va = pd.read_csv(ani_summary_path, sep="\t").rename(columns={
            "threshold_ani": "Umbral ANI Mash (%)", "n_pairs": "Pares", "n_with_fastani":
            "Con FastANI", "mean_difference": "Diferencia media (FastANI - Mash)",
            "median_abs_difference": "Diferencia absoluta mediana",
            "same_side_fraction": "Mismo lado del umbral"})  # fmt: skip
        blocks += [
            Block("h2", "Validacion de ANI: Mash frente a FastANI"),
            Block(
                "p",
                "Pares con ANI de Mash cerca de cada umbral, re-medidos con FastANI. "
                "'Mismo lado del umbral' es la fraccion de pares que ambos metodos ponen del "
                "mismo lado; detalle por par en derep/ani_validation.tsv.",
            ),
            Block("table", table=va),
        ]
        pairs_path = results / "derep/ani_validation.tsv"
        if pairs_path.exists():
            d = pd.read_csv(pairs_path, sep="\t")["difference"].dropna()
            if len(d):
                blocks.append(
                    Block(
                        "p",
                        f"En {int((d > 0).sum())} de {len(d)} pares FastANI dio un ANI mayor "
                        f"que Mash (diferencia mediana {d.median():+.2f} puntos). Los umbrales "
                        "de desreplicacion estan en la escala de Mash: si la diferencia es "
                        "sistematica, cada cluster es algo mas estricto que el mismo umbral "
                        "medido con FastANI. A esta resolucion ninguno de los dos metodos es "
                        "una referencia exacta.",
                    )
                )

    # --- Co-occurrence
    cooc_rows, cooc_tables = [], {}
    for path in sorted((results / "stats").glob("cooccurrence*.tsv")):
        name = path.stem.removeprefix("cooccurrence").removeprefix("_") or "all"
        c = pd.read_csv(path, sep="\t")
        cooc_tables[name] = c
        sig = c[c["q_value"] < alpha]
        cooc_rows.append(
            {
                "Conjunto": _set_label(name),
                "Genomas": int(totals.get(name, 0)),
                "Pares probados": len(c),
                f"q < {alpha:g}": len(sig),
                "Co-ocurrencia": int((sig["association"] == "co-occurrence").sum()),
                "Exclusion": int((sig["association"] == "exclusion").sum()),
                f"Ligados (<= {linkage_bp // 1000} kb)": int(sig["physically_linked"].sum()),
            }
        )
    order = {s: i for i, s in enumerate(sets)}
    cooc_summary = pd.DataFrame(cooc_rows)
    cooc_summary = cooc_summary.sort_values(
        "Conjunto", key=lambda s: s.map({_set_label(k): v for k, v in order.items()})
    )
    blocks += [
        Block("h2", "Co-ocurrencia"),
        Block(
            "p",
            "Fisher bilateral por par de tipos con 5-95 % de prevalencia, corregido por "
            "Benjamini-Hochberg. Un par esta 'ligado' si la distancia mediana entre los dos "
            "sistemas, en el mismo replicon, es menor o igual al umbral de ligamiento.",
        ),
        Block("table", table=cooc_summary),
    ]
    if cooc_set in cooc_tables:
        c = cooc_tables[cooc_set]
        free = c[(c["q_value"] < alpha) & ~c["physically_linked"]].sort_values("q_value")
        if len(free):
            t = free[["system_a", "system_b", "n_both", "odds_ratio", "q_value"]].copy()
            t["median_gap_bp"] = (free["median_gap_bp"] / 1000).round(0)
            t.columns = ["Tipo A", "Tipo B", "Genomas con ambos", "Odds ratio", "q",
                         "Distancia mediana (kb)"]  # fmt: skip
            blocks += [
                Block(
                    "p",
                    f"Pares asociados y no ligados en el conjunto {_set_label(cooc_set)}. "
                    "No implican interaccion funcional: los genomas siguen sin ser "
                    "independientes (no hay correccion filogenetica).",
                ),
                Block("table", table=t),
            ]

    # --- Figures and versions
    blocks.append(Block("h2", "Figuras"))
    captions = {
        "presence_heatmap": "Presencia/ausencia de tipos por genoma",
        "prevalence": "Prevalencia de los tipos mas frecuentes",
        "systems_per_genome": "Sistemas por genoma",
        "cooccurrence": "Co-ocurrencia",
    }
    for stem, caption in captions.items():
        if (results / f"figures/{stem}.png").exists():
            blocks.append(Block("img", f"figures/{stem}.png|{caption}"))
    blocks += [
        Block("h2", "Versiones"),
        Block("pre", (results / "versions.txt").read_text().strip()),
    ]
    return blocks


def to_markdown(blocks: list[Block]) -> str:
    """Render the report as Markdown (image paths relative to the run directory)."""
    out = []
    for b in blocks:
        if b.kind == "h1":
            out.append(f"# {b.text}")
        elif b.kind == "h2":
            out.append(f"## {b.text}")
        elif b.kind == "p":
            out.append(b.text)
        elif b.kind == "table":
            out.append(md_table(b.table))
        elif b.kind == "img":
            path, caption = b.text.split("|", 1)
            out.append(f"![{caption}]({path})")
        elif b.kind == "pre":
            out.append(f"```\n{b.text}\n```")
    return "\n\n".join(out) + "\n"


CSS = """
:root { --bg: #fcfcfb; --fg: #1f1f1d; --muted: #6b6a66; --rule: #e4e2dd; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #1b1b1a; --fg: #ececea; --muted: #a3a29e; --rule: #3a3936; }
}
body { background: var(--bg); color: var(--fg); font: 15px/1.55 system-ui, sans-serif;
       max-width: 960px; margin: 0 auto; padding: 24px 16px; }
h1 { font-size: 1.6em; } h2 { margin-top: 2em; border-bottom: 1px solid var(--rule); }
table { border-collapse: collapse; margin: 1em 0; font-size: 0.9em; display: block;
        overflow-x: auto; }
th, td { padding: 4px 10px; border-bottom: 1px solid var(--rule); text-align: left; }
th { color: var(--muted); font-weight: 600; }
figure { margin: 1.5em 0; } img { max-width: 100%; background: #fcfcfb; }
figcaption { color: var(--muted); font-size: 0.9em; }
pre { background: rgba(127,127,127,.1); padding: 12px; overflow-x: auto; }
"""


def to_html(blocks: list[Block], results: Path) -> str:
    """Render the report as one self-contained HTML file (figures embedded)."""
    body = []
    for b in blocks:
        if b.kind in ("h1", "h2"):
            body.append(f"<{b.kind}>{html.escape(b.text)}</{b.kind}>")
        elif b.kind == "p":
            body.append(f"<p>{html.escape(b.text)}</p>")
        elif b.kind == "table":
            body.append(b.table.map(_fmt).to_html(index=False, border=0, escape=True))
        elif b.kind == "img":
            path, caption = b.text.split("|", 1)
            data = base64.b64encode((results / path).read_bytes()).decode()
            body.append(
                f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}">'
                f"<figcaption>{html.escape(caption)}</figcaption></figure>"
            )
        elif b.kind == "pre":
            body.append(f"<pre>{html.escape(b.text)}</pre>")
    title = html.escape(blocks[0].text) if blocks else "pseudodefense"
    return (
        f'<!doctype html>\n<html lang="es"><head><meta charset="utf-8">'
        f'<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{title}</title><style>{CSS}</style></head><body>\n"
        + "\n".join(body)
        + "\n</body></html>\n"
    )
