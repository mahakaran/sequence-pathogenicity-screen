"""Annotate the saved UMAP by protein family and taxonomic domain.

Needs results/umap_coords.parquet from umap_projection.py. No network access.

Annotations (see annotate.py):
  * domain  - Bacteria / Eukaryota / Viruses / Archaea, from the NCBI phylum
  * family  - normalised protein name ("Superoxide dismutase [Cu-Zn] 2" ->
              "Superoxide dismutase [Cu-Zn]"). Proteins with uninformative names
              (uncharacterized, "Protein X") and the query ORFs get no family.

Colour is used one group at a time (small multiples): a scatter can only carry
~3 distinguishable hues, and there are ~77k families and 4 domains.

Outputs:
  results/umap_annotated_full.png   whole map: top families | domains
  results/umap_annotated_zoom.png   ORF neighbourhood: top families | domains
  results/umap_annotations.parquet  per-protein family/domain/coords
"""

import textwrap

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from annotate import DOMAINS, annotate
from umap_projection import BG_POINTS, HERE, INK, INK_2, INK_MUTED, ORF, PROJECT, SURFACE

HIGHLIGHT = "#2a78d6"   # palette slot 1; slot 2 (orange) marks the ORFs
N_FAMILIES = 8
LABEL_OFFSETS = {"ORF1": (-38, -14), "ORF2": (38, -6), "ORF3": (26, 22),
                 "ORF4": (38, 0), "ORF5": (-34, 14)}


def build_annotations():
    coords = pd.read_parquet(HERE / "results" / "umap_coords.parquet")
    n_ref = int((~coords.is_query).sum())
    meta = pq.read_table(PROJECT / "DS569k.parquet",
                         columns=["protein_name", "organism_name", "ncbi_taxonomy_phylum",
                                  "ncbi_taxonomy_class"])
    n_q = len(coords) - n_ref
    ann = annotate(meta.to_pandas())
    # Unnamed proteins and the ORFs get no family: they stay grey in family panels.
    coords["family"] = list(ann.family_named) + [None] * n_q
    coords["domain"] = list(ann.domain) + ["Query ORF"] * n_q
    coords.to_parquet(HERE / "results" / "umap_annotations.parquet", index=False)
    return coords


def style_axes(a):
    a.set_facecolor(SURFACE)
    a.set_xticks([]); a.set_yticks([])
    for sp in a.spines.values():
        sp.set_visible(False)


def facet(a, ref, q, mask, title, s_bg, s_hi, label_orfs):
    style_axes(a)
    a.scatter(ref.x, ref.y, s=s_bg, c=BG_POINTS, alpha=0.3, linewidths=0, rasterized=True)
    hi = ref[mask]
    a.scatter(hi.x, hi.y, s=s_hi, c=HIGHLIGHT, alpha=0.8, linewidths=0, rasterized=True)
    a.scatter(q.x, q.y, s=40, c=ORF, edgecolors=SURFACE, linewidths=1.5, zorder=5)
    if label_orfs:
        for _, r in q.iterrows():
            a.annotate(r.accession, (r.x, r.y), xytext=LABEL_OFFSETS.get(r.accession, (25, 15)),
                       textcoords="offset points", ha="center", va="center", fontsize=7.5,
                       fontweight="bold", color=INK, zorder=6,
                       arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.6, shrinkA=1, shrinkB=3),
                       bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.9))
    a.set_title(title, loc="left", fontsize=9, color=INK)


def figure(d, title, fname, s_bg, s_hi, label_orfs, window=None):
    ref, q = d[d.domain != "Query ORF"], d[d.domain == "Query ORF"]
    if window is not None:
        x0, x1, y0, y1 = window
        ref = ref[ref.x.between(x0, x1) & ref.y.between(y0, y1)]

    top = ref.family.dropna().value_counts().head(N_FAMILIES)

    fig = plt.figure(figsize=(20, 9.5), facecolor=SURFACE)
    outer = fig.add_gridspec(1, 2, width_ratios=[4, 2], wspace=0.06,
                             left=0.01, right=0.99, top=0.83, bottom=0.08)
    gf = outer[0].subgridspec(2, 4, hspace=0.18, wspace=0.04)
    gd = outer[1].subgridspec(2, 2, hspace=0.18, wspace=0.04)

    axes = []
    for k, (fam_name, n) in enumerate(top.items()):
        a = fig.add_subplot(gf[k // 4, k % 4])
        name = textwrap.shorten(fam_name, 70, placeholder="…")
        t = "\n".join(textwrap.wrap(name, 36)[:2]) + f"\nn={n:,}"
        facet(a, ref, q, (ref.family == fam_name).values, t, s_bg, s_hi, label_orfs and k == 0)
        axes.append(a)
    for k, dom_name in enumerate(DOMAINS):
        a = fig.add_subplot(gd[k // 2, k % 2])
        m = (ref.domain == dom_name).values
        facet(a, ref, q, m, f"{dom_name}\nn={m.sum():,} ({m.mean():.1%})", s_bg, s_hi,
              label_orfs and k == 0)
        axes.append(a)
    if window is not None:
        for a in axes:
            a.set_xlim(window[0], window[1]); a.set_ylim(window[2], window[3])

    fig.text(0.01, 0.96, title, fontsize=14, color=INK, ha="left")
    fig.text(0.01, 0.925, f"Protein family: {N_FAMILIES} largest in view (from normalised "
             "protein names; unnamed proteins left grey)", fontsize=10.5, color=INK_2)
    fig.text(0.675, 0.925, "Taxonomic domain", fontsize=10.5, color=INK_2)
    handles = [Line2D([], [], marker="o", ls="", ms=6, mfc=HIGHLIGHT, mec="none",
                      label="Proteins in the panel's group"),
               Line2D([], [], marker="o", ls="", ms=6, mfc=BG_POINTS, mec="none",
                      label="All other Swiss-Prot proteins"),
               Line2D([], [], marker="o", ls="", ms=7, mfc=ORF, mec=SURFACE, mew=1.5,
                      label="Query ORF")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=10,
               labelcolor=INK_2)
    out = HERE / "results" / fname
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out}")


def main():
    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})
    d = build_annotations()
    q = d[d.domain == "Query ORF"]
    ref = d[d.domain != "Query ORF"]
    print(f"Named family for {ref.family.notna().mean():.1%} of Swiss-Prot proteins")

    figure(d, "Swiss-Prot UMAP (ESM-2 8M + ProteinCLIP) — all 569,192 proteins",
           "umap_annotated_full.png", s_bg=0.03, s_hi=0.3, label_orfs=False)

    pad_x = max(np.ptp(q.x) * 0.25, 1.5); pad_y = max(np.ptp(q.y) * 0.25, 1.5)
    win = (q.x.min() - pad_x, q.x.max() + pad_x, q.y.min() - pad_y, q.y.max() + pad_y)
    n_in = int((ref.x.between(win[0], win[1]) & ref.y.between(win[2], win[3])).sum())
    figure(d, f"Zoom on the query ORFs — {n_in:,} proteins in view",
           "umap_annotated_zoom.png", s_bg=0.8, s_hi=3, label_orfs=True, window=win)


if __name__ == "__main__":
    main()
