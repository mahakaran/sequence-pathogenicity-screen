"""2D UMAP of all DS569k (Swiss-Prot) ProteinCLIP embeddings plus the query ORFs.

Step 1 (slow, once): fit UMAP on all 569k reference vectors + ORF vectors and
                     save coordinates to results/umap_coords.parquet.
Step 2 (fast):       plot. Re-run with --plot-only to restyle without refitting.

Usage:
    python umap_projection.py               # fit + plot
    python umap_projection.py --plot-only   # re-plot from saved coordinates
"""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
COORDS = HERE / "results" / "umap_coords.parquet"
FIG = HERE / "results" / "umap_swissprot_orfs.png"

# Reference palette (dataviz skill): scatter is capped at 3 categorical slots.
SURFACE, INK, INK_2, INK_MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
BG_POINTS = "#b9b8b3"
SOD = "#2a78d6"       # slot 1: Cu-Zn superoxide dismutases
ORF = "#eb6834"       # slot 2: query ORFs


def embed_orfs(query):
    """Embed ORFs in a child process.

    torch and faiss each ship their own OpenMP runtime on macOS; loading both
    into one process segfaults. So torch only ever runs in the subprocess.
    """
    import subprocess, sys
    recs = list(SeqIO.parse(query, "fasta"))
    seqs = [str(r.seq).rstrip("*").upper() for r in recs]
    names = [r.id.split("|")[-1].split(":")[0] for r in recs]   # "ORF1", ...
    out = HERE / "results" / "orf_proteinclip.npy"
    out.parent.mkdir(parents=True, exist_ok=True)
    code = ("import sys, numpy as np; sys.path.insert(0, sys.argv[1]);"
            "from esm_search import Embedder, ProteinCLIP;"
            "seqs = sys.argv[3:];"
            "np.save(sys.argv[2], ProteinCLIP()(Embedder().embed(seqs)))")
    subprocess.run([sys.executable, "-c", code, str(HERE), str(out), *seqs],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return names, seqs, np.load(out)


def fit(args):
    import umap
    t = pq.read_table(args.refs, columns=["accession", "protein_name", "organism_name",
                                          "ncbi_taxonomy_phylum", "sequence_length",
                                          "embedding"])
    meta = t.drop(["embedding"]).to_pandas()
    ref = np.asarray(t.column("embedding").combine_chunks().values,
                     dtype=np.float32).reshape(len(meta), -1)
    ref = ref / np.linalg.norm(ref, axis=1, keepdims=True)   # arrow buffer is read-only

    names, seqs, q = embed_orfs(args.query)
    X = np.ascontiguousarray(np.vstack([ref, q.astype(np.float32)]))

    # kNN graph with FAISS HNSW (multi-threaded, ~1 min) instead of UMAP's
    # built-in NN-descent, which is the slow, mostly single-threaded step here.
    import faiss
    K = 15
    t0 = time.time()
    index = faiss.IndexHNSWFlat(X.shape[1], 32, faiss.METRIC_INNER_PRODUCT)
    index.hnsw.efConstruction = 80
    index.add(X)
    index.hnsw.efSearch = 64
    sim, knn_idx = index.search(X, K)
    # Unit vectors: euclidean distance = sqrt(2 - 2*cos). UMAP expects self as 1st neighbour.
    knn_dist = np.sqrt(np.clip(2 - 2 * sim, 0, None)).astype(np.float32)
    # Identical vectors (duplicate sequences) can outrank self; force self into column 0.
    rows = np.arange(len(X))
    not_self = knn_idx[:, 0] != rows
    knn_idx[not_self, 0], knn_dist[not_self, 0] = rows[not_self], 0.0
    print(f"kNN graph in {time.time() - t0:.0f}s ({not_self.mean():.1%} rows had a "
          f"duplicate ranked above self)", flush=True)

    print(f"Fitting UMAP on {X.shape[0]:,} x {X.shape[1]} ...", flush=True)
    t0 = time.time()
    reducer = umap.UMAP(n_neighbors=K, min_dist=0.1, metric="euclidean",
                        precomputed_knn=(knn_idx, knn_dist, None),
                        # spectral init (default) is a single-threaded eigensolve
                        # over the whole 569k graph and stalls for 30+ min
                        init="pca",
                        low_memory=True, n_jobs=-1, verbose=True)
    Y = reducer.fit_transform(X)
    print(f"UMAP done in {(time.time() - t0) / 60:.1f} min", flush=True)

    orf_meta = pd.DataFrame({"accession": names, "protein_name": "query ORF",
                             "organism_name": None, "ncbi_taxonomy_phylum": None,
                             "sequence_length": [len(s) for s in seqs]})
    out = pd.concat([meta, orf_meta], ignore_index=True)
    out["is_query"] = [False] * len(meta) + [True] * len(names)
    out["x"], out["y"] = Y[:, 0], Y[:, 1]
    COORDS.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(COORDS, index=False)
    print(f"Wrote {COORDS}")


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = pd.read_parquet(COORDS)
    ref, q = d[~d.is_query], d[d.is_query]
    sod = ref.protein_name.str.contains(r"superoxide dismutase \[Cu-Zn\]", case=False, na=False)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "text.color": INK, "axes.labelcolor": INK_2})
    fig, (ax, axz) = plt.subplots(1, 2, figsize=(16, 8), facecolor=SURFACE,
                                  gridspec_kw={"width_ratios": [1, 1]})

    def base(a, sub, s_bg, s_sod):
        a.set_facecolor(SURFACE)
        a.scatter(sub.x, sub.y, s=s_bg, c=BG_POINTS, alpha=0.35, linewidths=0,
                  rasterized=True)
        ss = sub[sod.reindex(sub.index, fill_value=False)]
        a.scatter(ss.x, ss.y, s=s_sod, c=SOD, alpha=0.9, linewidths=0, rasterized=True)
        a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)

    def orfs(a, labels=True):
        a.scatter(q.x, q.y, s=90, c=ORF, edgecolors=SURFACE, linewidths=2, zorder=5)
        if not labels:
            return
        # Fixed label offsets (points), with leader lines; ORF1/ORF2 nearly coincide
        offsets = {"ORF1": (-70, -25), "ORF2": (70, -10), "ORF3": (40, 40),
                   "ORF4": (70, 0), "ORF5": (-60, 25)}
        for _, r in q.iterrows():
            off = offsets.get(r.accession, (50, 30))
            a.annotate(r.accession, (r.x, r.y), xytext=off, textcoords="offset points",
                       ha="center", va="center", fontsize=10, fontweight="bold", color=INK,
                       zorder=6, arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8,
                                                 shrinkA=2, shrinkB=5),
                       bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none", alpha=0.9))

    # Left: whole map
    base(ax, ref, 0.05, 1.5)
    orfs(ax, labels=False)
    ax.annotate("query ORFs", (q.x.mean(), q.y.max()), xytext=(0, 40),
                textcoords="offset points", ha="center", fontsize=10, fontweight="bold",
                color=INK, arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8))
    ax.set_title(f"All {len(ref):,} Swiss-Prot proteins + 5 query ORFs",
                 loc="left", fontsize=12, color=INK)

    # Right: zoom around the ORFs (padded bounding box)
    pad_x = max(np.ptp(q.x) * 0.25, 1.5); pad_y = max(np.ptp(q.y) * 0.25, 1.5)
    x0, x1 = q.x.min() - pad_x, q.x.max() + pad_x
    y0, y1 = q.y.min() - pad_y, q.y.max() + pad_y
    win = ref[(ref.x.between(x0, x1)) & (ref.y.between(y0, y1))]
    base(axz, win, 1.2, 6)
    orfs(axz)
    axz.set_xlim(x0, x1); axz.set_ylim(y0, y1)
    axz.set_title(f"Zoom on the ORFs ({len(win):,} proteins in view)",
                  loc="left", fontsize=12, color=INK)
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               ec=INK_MUTED, lw=1, ls="--"))

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", ls="", ms=6, mfc=BG_POINTS, mec="none", label="Swiss-Prot protein"),
        Line2D([], [], marker="o", ls="", ms=6, mfc=SOD, mec="none",
               label=f"Cu-Zn superoxide dismutase (n={int(sod.sum()):,})"),
        Line2D([], [], marker="o", ls="", ms=9, mfc=ORF, mec=SURFACE, mew=2, label="Query ORF"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=10, labelcolor=INK_2)
    fig.suptitle("UMAP of ESM-2 (8M) + ProteinCLIP embeddings", x=0.02, ha="left",
                 fontsize=14, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(FIG, dpi=200, facecolor=SURFACE)
    print(f"Wrote {FIG}")

    # Nearest-in-2D check: what surrounds each ORF on the map
    for _, r in q.iterrows():
        dist = np.hypot(ref.x - r.x, ref.y - r.y)
        nn = ref.loc[dist.nsmallest(20).index, "protein_name"].value_counts().head(3)
        print(f"{r.accession}: 20 nearest on map -> " + "; ".join(f"{n}× {k[:45]}" for k, n in nn.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", type=Path, default=PROJECT / "screening-exercise-orfs.fasta")
    ap.add_argument("--refs", type=Path, default=PROJECT / "DS569k.parquet")
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()
    if not args.plot_only:
        fit(args)
    plot()


if __name__ == "__main__":
    main()
