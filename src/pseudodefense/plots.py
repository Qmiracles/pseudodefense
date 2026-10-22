"""Static figures (PNG + SVG) for the README and the report.

Colors follow one rule per job: a single blue hue for presence and magnitude,
a blue <-> red diverging scale with a neutral gray midpoint for log2 odds ratios.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
ABSENT = "#f0efec"
BLUE = "#2a78d6"
BLUE_DARK = "#184f95"
RED = "#e34948"
NEUTRAL = "#f0efec"

STYLE = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK_SECONDARY,
    "axes.edgecolor": AXIS,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "svg.fonttype": "none",
}

DIVERGING = LinearSegmentedColormap.from_list("blue_gray_red", [RED, NEUTRAL, BLUE])


def _save(fig: plt.Figure, out: Path) -> None:
    """Write ``out`` (PNG) and an SVG next to it."""
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def presence_heatmap(matrix: pd.DataFrame, out: Path, title: str) -> None:
    """Clustered genome x system heatmap (Jaccard distance, average linkage)."""
    data = matrix.loc[:, matrix.sum(axis=0) > 0]
    with plt.rc_context(STYLE):
        width = min(4 + 0.18 * data.shape[1], 30)
        height = min(3 + 0.16 * data.shape[0], 30)
        grid = sns.clustermap(
            data.astype(bool),
            metric="jaccard",
            method="average",
            cmap=ListedColormap([ABSENT, BLUE]),
            linewidths=0.4 if data.shape[0] <= 100 else 0,
            linecolor=SURFACE,
            xticklabels=True,
            yticklabels=data.shape[0] <= 100,
            figsize=(width, height),
            dendrogram_ratio=(0.08, 0.08),
            cbar_pos=None,
            tree_kws={"colors": MUTED, "linewidths": 0.6},
        )
        grid.ax_heatmap.set_xlabel("Defense system (DefenseFinder type)")
        grid.ax_heatmap.set_ylabel("Genome")
        grid.ax_heatmap.tick_params(axis="x", labelsize=7, length=0)
        grid.ax_heatmap.tick_params(axis="y", labelsize=6, length=0)
        # Title and legend sit above the column dendrogram, not on top of it
        grid.figure.suptitle(title, x=0.02, y=1.06, ha="left", fontweight="bold", fontsize=11)
        handles = [
            plt.Rectangle((0, 0), 1, 1, color=BLUE, label="Detected"),
            plt.Rectangle((0, 0), 1, 1, color=ABSENT, label="Not detected"),
        ]
        grid.figure.legend(
            handles=handles,
            loc="lower right",
            bbox_to_anchor=(0.98, 1.01),
            ncol=2,
            frameon=False,
            fontsize=8,
        )
        _save(grid.figure, out)


def prevalence_bars(prev: pd.DataFrame, out: Path, title: str, top_n: int = 30) -> None:
    """Horizontal bars of prevalence with 95 % Wilson intervals (top ``top_n``)."""
    data = prev.head(top_n).iloc[::-1]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6.5, 0.6 + 0.24 * len(data)))
        y = np.arange(len(data))
        ax.barh(y, data["prevalence"] * 100, height=0.6, color=BLUE, zorder=2)
        ax.errorbar(
            data["prevalence"] * 100,
            y,
            xerr=[
                (data["prevalence"] - data["ci95_low"]) * 100,
                (data["ci95_high"] - data["prevalence"]) * 100,
            ],
            fmt="none",
            ecolor=BLUE_DARK,
            elinewidth=1,
            capsize=0,
            zorder=3,
        )
        ax.set_yticks(y, data["system"])
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(0, 100)
        ax.set_xlabel(f"Genomes carrying the system (%, n = {int(data['total_genomes'].iloc[0])})")
        ax.xaxis.grid(True, zorder=0)
        ax.spines["left"].set_visible(False)
        ax.set_title(title)
        ax.text(
            0,
            -0.6 / len(data) - 0.06,
            "Bars: % of genomes. Whiskers: 95 % Wilson interval.",
            transform=ax.transAxes,
            color=MUTED,
            fontsize=7,
            va="top",
        )
        _save(fig, out)


def systems_histogram(counts: pd.DataFrame, out: Path, title: str) -> None:
    """Distribution of the number of defense systems per genome."""
    values = counts["n_systems"]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6, 3.2))
        bins = np.arange(values.min() - 0.5, values.max() + 1.5, 1)
        ax.hist(values, bins=bins, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=2)
        ax.yaxis.get_major_locator().set_params(integer=True)
        ax.xaxis.get_major_locator().set_params(integer=True)
        top = ax.get_ylim()[1] * 1.15  # headroom so the median label clears the bars
        ax.set_ylim(0, top)
        ax.axvline(values.median(), color=INK_SECONDARY, linewidth=1, linestyle="--", zorder=3)
        ax.text(
            values.median(),
            top,
            f" median {values.median():g}",
            color=INK_SECONDARY,
            va="top",
            fontsize=8,
        )
        ax.set_xlabel("Defense systems per genome")
        ax.set_ylabel("Genomes")
        ax.yaxis.grid(True, zorder=0)
        ax.set_title(title)
        _save(fig, out)


def cooccurrence_heatmap(
    cooc: pd.DataFrame, out: Path, title: str, alpha: float = 0.05, max_systems: int = 40
) -> None:
    """Symmetric matrix of Haldane-corrected log2 odds ratios; dots mark q < alpha.

    Only the ``max_systems`` systems involved in the most pairs with the lowest
    q-values are shown, to keep labels legible.
    """
    if cooc.empty:
        with plt.rc_context(STYLE):
            fig, ax = plt.subplots(figsize=(5, 2))
            ax.axis("off")
            ax.text(0.5, 0.5, "No system pairs met the prevalence filter", ha="center", color=MUTED)
            ax.set_title(title)
            _save(fig, out)
        return
    ranked = pd.concat(
        [
            cooc[["system_a", "q_value"]].rename(columns={"system_a": "s"}),
            cooc[["system_b", "q_value"]].rename(columns={"system_b": "s"}),
        ]
    )
    keep = ranked.groupby("s")["q_value"].min().sort_values().head(max_systems).index
    keep = sorted(keep)
    lor = pd.DataFrame(np.nan, index=keep, columns=keep)
    sig = pd.DataFrame(False, index=keep, columns=keep)
    linked = pd.DataFrame(False, index=keep, columns=keep)
    has_linkage = "physically_linked" in cooc.columns
    for r in cooc.itertuples(index=False):
        if r.system_a in lor.index and r.system_b in lor.index:
            for i, j in ((r.system_a, r.system_b), (r.system_b, r.system_a)):
                lor.loc[i, j] = r.log2_or_haldane
                sig.loc[i, j] = r.q_value < alpha
                linked.loc[i, j] = has_linkage and bool(r.physically_linked)
    limit = max(1.0, float(np.nanmax(np.abs(lor.to_numpy()))))
    # The matrix is symmetric: show the lower triangle only
    upper = np.triu(np.ones(lor.shape, dtype=bool))
    lower_sig = sig.to_numpy() & ~upper & ~linked.to_numpy()
    lower_sig_linked = sig.to_numpy() & ~upper & linked.to_numpy()
    is_sig = cooc["q_value"] < alpha
    n_sig = int(is_sig.sum())
    n_linked = int((is_sig & cooc["physically_linked"]).sum()) if has_linkage else 0
    with plt.rc_context(STYLE):
        size = 2.5 + 0.22 * len(keep)
        fig, ax = plt.subplots(figsize=(size + 1.2, size))
        sns.heatmap(
            lor,
            mask=upper,
            cmap=DIVERGING,
            vmin=-limit,
            vmax=limit,
            center=0,
            square=True,
            linewidths=0.4,
            linecolor=SURFACE,
            cbar_kws={"label": "log2 odds ratio (Haldane-corrected)", "shrink": 0.6},
            ax=ax,
        )
        ys, xs = np.nonzero(lower_sig)
        ax.scatter(xs + 0.5, ys + 0.5, s=10, color=INK, zorder=3)
        ys, xs = np.nonzero(lower_sig_linked)
        ax.scatter(
            xs + 0.5, ys + 0.5, s=28, facecolors="none", edgecolors=INK, linewidths=1, zorder=3
        )
        ax.tick_params(length=0, labelsize=7)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(title)
        caption = (
            "Blue: found together more often than expected. Red: less often.\n"
            f"Filled dot: q < {alpha} (Fisher, Benjamini-Hochberg). "
            f"Open ring: q < {alpha} but the two systems sit in the same locus "
            "(physically linked), so co-occurrence is expected.\n"
            f"{len(cooc)} pairs tested, {n_sig} with q < {alpha}, "
            f"{n_linked} of them physically linked."
        )
        fig.text(0.01, -0.02, caption, color=MUTED, fontsize=7, va="top")
        _save(fig, out)
