"""Plot ESM-2 650M pseudo-log-likelihood: ORFs against length-matched controls.

Panel A  per-ORF strip plot - each ORF against its own matched controls
Panel B  all controls pooled, as distributions, with the ORFs marked
Panel C  PLL vs length, showing the length trend the matching controls for

Usage:  python pll_figure.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
FIG = RES / "pll_controls.png"

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
# validated 4-colour set (same as the UMAP domains), ORF in the figure accent
COLORS = {"uniprot": "#2a78d6", "shuffled": "#eda100",
          "comp_match": "#008300", "uniform": "#e34948"}
ORF_COLOR = "#eb6834"
NICE = {"uniprot": "Real Swiss-Prot (length-matched)",
        "shuffled": "Those proteins, shuffled",
        "comp_match": "ORF residues shuffled (same composition)",
        "uniform": "Uniform random over 20 aa"}
ORDER = ["uniprot", "shuffled", "comp_match", "uniform"]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    d = pd.read_parquet(RES / "pll_scores.parquet")
    orfs = d[d.group == "ORF"].set_index("orf")
    ctl = d[d.group != "ORF"]
    names = list(orfs.index)

    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})
    fig = plt.figure(figsize=(16, 12.4), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1],
                          top=0.852, bottom=0.085, left=0.065, right=0.985,
                          hspace=0.30, wspace=0.20)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    def style(ax):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_2, labelsize=9)
        ax.grid(axis="y", color=GRID, lw=0.7, alpha=0.8)
        ax.set_axisbelow(True)

    # ---- A: per ORF -------------------------------------------------------
    style(axA)
    rng = np.random.default_rng(0)
    width = 0.19
    for i, orf in enumerate(names):
        for j, g in enumerate(ORDER):
            v = ctl[(ctl.orf == orf) & (ctl.group == g)].pll.values
            if not len(v):
                continue
            x = i + (j - 1.5) * width
            axA.scatter(x + rng.normal(0, 0.022, len(v)), v, s=26, color=COLORS[g],
                        alpha=0.75, linewidths=0, zorder=3)
            axA.plot([x - width * 0.38, x + width * 0.38], [v.mean()] * 2,
                     color=COLORS[g], lw=2.2, zorder=4)
        y = orfs.loc[orf, "pll"]
        axA.scatter([i], [y], s=210, marker="*", color=ORF_COLOR,
                    edgecolors="white", linewidths=1.3, zorder=6)
        # label up and left of the star: the comp_match group sits immediately right
        axA.annotate(f"{y:+.2f}", (i - 1.5 * width, y), xytext=(-6, 7),
                     textcoords="offset points", fontsize=9, color=INK,
                     fontweight="bold", ha="right", zorder=8,
                     bbox=dict(boxstyle="round,pad=0.16", fc=SURFACE, ec="none", alpha=0.9))
    axA.set_xticks(range(len(names)))
    axA.set_xticklabels([f"{o}\n{int(orfs.loc[o, 'length'])} aa" for o in names], fontsize=10)
    axA.set_ylabel("pseudo-log-likelihood per residue", fontsize=10.5, color=INK_2)
    axA.set_title("A · Each ORF against its own length-matched controls",
                  loc="left", fontsize=13, color=INK, pad=34)
    # a little headroom for the legend, without leaving a band of empty axis
    lo, hi = axA.get_ylim()
    axA.set_ylim(lo - 0.04 * (hi - lo), hi + 0.08 * (hi - lo))
    h = [Line2D([], [], marker="o", ls="", ms=8, mfc=COLORS[g], mec="none", label=NICE[g])
         for g in ORDER]
    h.append(Line2D([], [], marker="*", ls="", ms=14, mfc=ORF_COLOR, mec="white",
                    label="query ORF"))
    axA.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 1.085),
               frameon=False, fontsize=9.2, labelcolor=INK_2, ncol=5,
               handletextpad=0.5, columnspacing=1.6)

    # ---- B: pooled distributions -----------------------------------------
    style(axB)
    for g in ORDER:
        v = ctl[ctl.group == g].pll.values
        if len(v):
            axB.hist(v, bins=26, alpha=0.55, color=COLORS[g], label=NICE[g],
                     edgecolor="none")
    # ORF markers: ticks below the axis with staggered labels, so five lines within
    # ~0.35 PLL units do not collide with each other or with the bars.
    top = axB.get_ylim()[1]
    axB.set_ylim(0, top * 1.40)
    for k, orf in enumerate(sorted(names, key=lambda o: orfs.loc[o, "pll"])):
        xv = orfs.loc[orf, "pll"]
        tier = top * (1.03 + 0.11 * (k % 3))
        # draw to the label's own tier so the rules do not run past the staggered boxes
        axB.vlines(xv, 0, tier, color=ORF_COLOR, lw=1.3, alpha=0.85, zorder=5)
        axB.annotate(orf, (xv, tier), fontsize=8, color=INK, ha="center", va="bottom",
                     zorder=7,
                     bbox=dict(boxstyle="round,pad=0.17", fc=SURFACE, ec=GRID, lw=0.5))
    axB.set_xlabel("pseudo-log-likelihood per residue", fontsize=10.5, color=INK_2)
    axB.set_ylabel("sequences", fontsize=10.5, color=INK_2)
    axB.set_title("B · All controls pooled", loc="left", fontsize=13, color=INK, pad=10)
    axB.legend(frameon=True, facecolor=SURFACE, edgecolor=GRID, fontsize=7.6,
               labelcolor=INK_2, framealpha=0.95, loc="center left")

    # ---- C: PLL vs length -------------------------------------------------
    style(axC)
    for g in ORDER:
        s = ctl[ctl.group == g]
        axC.scatter(s.length, s.pll, s=22, color=COLORS[g], alpha=0.6,
                    linewidths=0, label=NICE[g])
    axC.scatter(orfs.length, orfs.pll, s=200, marker="*", color=ORF_COLOR,
                edgecolors="white", linewidths=1.3, zorder=6)
    # ORF1 and ORF4 are both 96 aa and nearly equal PLL, so fixed offsets separate them
    OFF = {"ORF1": (12, 9), "ORF4": (12, -14), "ORF2": (-14, 9),
           "ORF3": (12, 4), "ORF5": (12, 4)}
    for orf in names:
        axC.annotate(orf, (orfs.loc[orf, "length"], orfs.loc[orf, "pll"]),
                     xytext=OFF.get(orf, (10, 4)), textcoords="offset points",
                     fontsize=8.4, color=INK, fontweight="bold", zorder=7,
                     bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none",
                               alpha=0.85))
    axC.set_xlabel("sequence length (aa)", fontsize=10.5, color=INK_2)
    axC.set_ylabel("pseudo-log-likelihood per residue", fontsize=10.5, color=INK_2)
    axC.set_title("C · PLL against length", loc="left", fontsize=13, color=INK, pad=10)

    # ---- percentile summary ----------------------------------------------
    # percentile footer, wrapped onto two lines so it stays inside the canvas
    lines = []
    for orf in names:
        y = orfs.loc[orf, "pll"]
        bits = []
        for g in ORDER:
            v = ctl[(ctl.orf == orf) & (ctl.group == g)].pll.values
            if len(v):
                bits.append(f"{g} {100 * (v < y).mean():.0f}%")
        lines.append(f"{orf}: " + ", ".join(bits))
    fig.text(0.065, 0.042,
             "Percentile of each ORF within its own control groups "
             "(fraction of controls scoring lower)",
             fontsize=8.6, color=INK_2, ha="left", fontweight="bold")
    half = int(np.ceil(len(lines) / 2))
    for r, chunk in enumerate((lines[:half], lines[half:])):
        fig.text(0.065, 0.024 - r * 0.016, "     |     ".join(chunk),
                 fontsize=8.2, color=INK_MUTED, ha="left")

    n = ctl.groupby("group").size().min()
    fig.suptitle("Are the ORFs protein-like? ESM-2 650M pseudo-log-likelihood",
                 x=0.065, y=0.978, ha="left", fontsize=17, color=INK)
    fig.text(0.065, 0.950,
             f"Masked-position PLL, mean per residue. {n} controls per group per ORF, each "
             "length-matched to the ORF it is compared against. Higher = more protein-like "
             "to the model.",
             ha="left", fontsize=10.3, color=INK_2, va="top")
    fig.savefig(FIG, dpi=200, facecolor=SURFACE)
    print(f"Wrote {FIG}")

    print("\nGroup means (pooled):")
    print(ctl.groupby("group").pll.agg(["count", "mean", "std"]).round(3))
    print("\nPer-ORF percentile within each control group:")
    for ln in lines:
        print("  " + ln)


if __name__ == "__main__":
    main()
