"""
Two figures over the Leiden-start sweeps of Experiment 1 (all pairs, trees and
mappers), one panel per reference, a box per criterion (set) and the sweeps as
dots:

    selected_k_per_reference     the K each criterion set selects (selected_k.csv,
                                 01_select.py), the annotation's T dashed
    score_range_per_reference    each criterion score's range over K
                                 (max - min per sweep, all_k_per_run.csv.gz)

Both are written as PNG and PDF to ABLATION_CRITERIA.
"""

from __future__ import annotations
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from spatialaim_experiments import paths, style

style.use_paper_rcparams()

# The sets in display order with their tick labels and colours: each criterion
# in its own colour, the pairs grey, all three in ink.
SETS = {
    "recon": ("R", style.C_RECON),
    "spatial": ("S", style.C_SPATIAL),
    "coherence": ("C", style.C_COHERENCE),
    "recon+spatial": ("R+S", style.GREY),
    "recon+coherence": ("R+C", style.GREY),
    "spatial+coherence": ("S+C", style.GREY),
    "all": ("all", style.INK),
}
CRITERIA = ["recon", "spatial", "coherence"]
K_TICKS = [1, 2, 4, 8, 16, 32, 64]


def leiden_sweeps() -> pd.DataFrame:
    """Every K row of the Leiden-start sweeps, with the three criterion scores."""
    allk = pd.read_csv(paths.EXP1_SUMMARY / "all_k_per_run.csv.gz",
                       usecols=["k", "run_id", "pair_id", "sc", "start_tag", *CRITERIA])
    return allk[allk["start_tag"] == "leiden"]


def annotation_types() -> pd.Series:
    """Per reference the sorted numbers of types T of its annotations: the first
    K of each annotation-start sweep."""
    allk = pd.read_csv(paths.EXP1_SUMMARY / "all_k_per_run.csv.gz",
                       usecols=["k", "sc", "start_tag"])
    allk = allk[allk["start_tag"] != "leiden"]
    return allk.groupby(["sc", "start_tag"])["k"].max().groupby("sc").unique().map(sorted)


def score_ranges(allk: pd.DataFrame) -> pd.DataFrame:
    """Per sweep and criterion the score's max - min over K, long, as ``value``."""
    by_run = allk.groupby(["run_id", "pair_id", "sc"])[CRITERIA]
    wide = (by_run.max() - by_run.min()).reset_index()
    return wide.melt(id_vars=["run_id", "pair_id", "sc"], var_name="criteria", value_name="value")


def panel(ax, df: pd.DataFrame, sets: list[str], rng: np.random.Generator) -> None:
    """One reference: per set a box and the dots of ``df["value"]``."""
    for x, criteria in enumerate(sets):
        colour = SETS[criteria][1]
        v = df.loc[df["criteria"] == criteria, "value"].dropna().to_numpy(float)
        ax.scatter(x + rng.uniform(-0.28, 0.28, len(v)), v, s=2.5, color=colour,
                   alpha=0.25, linewidths=0, zorder=1)
        ax.boxplot(v, positions=[x], widths=0.6, showfliers=False, zorder=2,
                   medianprops=dict(color=colour, lw=1.4),
                   boxprops=dict(color=colour, lw=0.8), whiskerprops=dict(color=colour, lw=0.8),
                   capprops=dict(color=colour, lw=0.8))
    ax.set_xticks(range(len(sets)), [SETS[s][0] for s in sets],
                  rotation=90 if len(sets) > 3 else 0)
    style.bare(ax)


def figure(values: pd.DataFrame, sets: list[str], ylabel: str, note: str, stem: str,
           decorate) -> None:
    """One panel per reference of ``values`` (``sc``, ``criteria``, ``value``);
    ``decorate(ax, reference)`` adds the axis scale and reference lines."""
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(2, 5, figsize=(7.2, 3.6), sharey=True)
    for ax, (reference, label) in zip(axes.flat, style.REFERENCE_LABELS.items()):
        df = values[values["sc"] == reference]
        panel(ax, df, sets, rng)
        n_pairs = df["pair_id"].nunique()
        ax.set_title(f"{label}\n{n_pairs} pair{'s' if n_pairs > 1 else ''}, "
                     f"{df['run_id'].nunique()} sweeps", fontsize=7.5)
        decorate(ax, reference)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    fig.text(0.99, 0.005, note, ha="right", va="bottom", fontsize=7, color=style.GREY)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(paths.ABLATION_CRITERIA / f"{stem}.{ext}",
                    dpi=style.PREVIEW_DPI, facecolor="white")
    plt.close(fig)


def main() -> None:
    selected = pd.read_csv(paths.ABLATION_CRITERIA / "selected_k.csv").rename(columns={"k": "value"})
    types = annotation_types()

    def k_axis(ax, reference):
        ax.set_yscale("log", base=2)
        ax.set_yticks(K_TICKS, [str(k) for k in K_TICKS])
        ax.minorticks_off()
        ax.set_ylim(0.8, 100)
        for t in types.get(reference, []):
            ax.axhline(t, color=style.GREY, lw=0.7, ls=(0, (3, 2)), zorder=0)

    figure(selected, list(SETS), "selected K",
           "R recon, S spatial, C coherence; dashed: annotation's T",
           "selected_k_per_reference", k_axis)

    ranges = score_ranges(leiden_sweeps())
    top = 1.05 * ranges["value"].max()
    figure(ranges, CRITERIA, "score range over K",
           "range = max − min of the score over a sweep's K; R recon, S spatial, C coherence",
           "score_range_per_reference", lambda ax, reference: ax.set_ylim(0, top))


if __name__ == "__main__":
    main()
