"""UMAP of all Swiss-Prot ESM-2 650M embeddings + the query ORFs, as a 4-panel figure.

Stages (each cached, so re-runs are cheap):
  embed   embed the 5 ORFs with ESM-2 650M using the pooling convention that
          matches downloaded_embeddings/embeddings_esm2_650m.h5 (verified first)
  umap    PCA -> FAISS HNSW kNN -> UMAP over 540,881 proteins + 5 ORFs
  annot   per-protein taxonomic domain and InterPro homologous superfamily
  plot    one figure, exactly four subplots

Usage:
    python uniprot_650m_map.py verify     # which pooling reproduces the h5?
    python uniprot_650m_map.py embed
    python uniprot_650m_map.py umap
    python uniprot_650m_map.py plot
    python uniprot_650m_map.py all
"""

import csv
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
H5 = HERE / "downloaded_embeddings" / "embeddings_esm2_650m.h5"
MODEL_DIR = HERE / "model_650m"
XREFS = HERE / "annotations" / "sprot_xrefs.tsv.gz"
ENTRY_LIST = HERE / "annotations" / "entry.list"
RES = HERE / "results"
ORF_EMB = RES / "orf_esm2_650m.npy"
MISSING_EMB = RES / "missing_650m.npy"
MISSING_ACC = RES / "missing_650m_accessions.txt"
ORF_IDS = RES / "orf_ids.txt"
COORDS = RES / "umap650_coords.parquet"
FIG = RES / "uniprot_650m_umap_4panel.png"
MAX_RES = 1022

# ---- colours -------------------------------------------------------------
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
UNASSIGNED = "#d8d7d2"
ORF_COLOR = "#eb6834"
# validated 4-colour set for the taxonomic domains (all-pairs, light surface)
DOMAIN_COLORS = {"Bacteria": "#2a78d6", "Eukaryota": "#eda100",
                 "Archaea": "#008300", "Viruses": "#e34948"}
N_SUPERFAMILIES = 20          # coloured from tab20; the rest fold into "Other"
N_ZOOM_SUPERFAMILIES = 20     # panel D ranks these within the zoom window, not globally


# =============================================================== embedding
def load_model():
    import torch
    from transformers import AutoModel, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(str(MODEL_DIR))
    model = AutoModel.from_pretrained(str(MODEL_DIR)).eval()
    dev = "mps" if torch.backends.mps.is_available() else "cpu"
    return tok, model.to(dev), dev


def pool_variants(tok, model, dev, seqs):
    """Return {variant_name: (n, 1280) array} for the plausible conventions."""
    import torch
    out = {k: [] for k in ("mean_nospecial", "mean_all", "cls")}
    with torch.no_grad():
        for s in seqs:
            enc = tok([s[:MAX_RES]], return_tensors="pt").to(dev)
            h = model(**enc).last_hidden_state[0]
            L = min(len(s), MAX_RES)
            out["mean_nospecial"].append(h[1:1 + L].mean(0).float().cpu().numpy())
            out["mean_all"].append(h.mean(0).float().cpu().numpy())
            out["cls"].append(h[0].float().cpu().numpy())
    return {k: np.vstack(v) for k, v in out.items()}


def cmd_verify(n_test=6):
    """Find which pooling convention reproduces the downloaded embeddings."""
    import h5py
    import pyarrow.parquet as pq
    t = pq.read_table(PROJECT / "DS569k.parquet", columns=["accession", "sequence"]).to_pandas()
    acc2seq = dict(zip(t.accession, t.sequence))
    with h5py.File(H5, "r") as f:
        cand = [a for a in t.accession if a in f and 40 <= len(acc2seq[a]) <= 250][:n_test]
        ref = np.vstack([f[a][()] for a in cand])
    tok, model, dev = load_model()
    print(f"device: {dev}; testing {len(cand)} proteins")
    best, best_cos = None, -1
    for name, V in pool_variants(tok, model, dev, [acc2seq[a] for a in cand]).items():
        cos = (V * ref).sum(1) / np.linalg.norm(V, axis=1) / np.linalg.norm(ref, axis=1)
        print(f"  {name:<16} cosine vs stored: min={cos.min():.6f} mean={cos.mean():.6f}")
        if cos.mean() > best_cos:
            best, best_cos = name, cos.mean()
    print(f"\n-> best match: {best} (mean cosine {best_cos:.6f})")
    if best_cos < 0.99:
        print("   WARNING: no convention reproduces the file closely; "
              "ORF vectors may not be directly comparable.")
    return best


def cmd_embed(convention=None):
    from Bio import SeqIO
    convention = convention or cmd_verify()
    recs = list(SeqIO.parse(PROJECT / "screening-exercise-orfs.fasta", "fasta"))
    names = [r.id.split("|")[-1].split(":")[0] for r in recs]
    seqs = [str(r.seq).rstrip("*").upper() for r in recs]
    tok, model, dev = load_model()
    V = pool_variants(tok, model, dev, seqs)[convention]
    RES.mkdir(exist_ok=True)
    np.save(ORF_EMB, V)
    ORF_IDS.write_text("\n".join(names))
    print(f"Embedded {len(names)} ORFs with '{convention}' -> {ORF_EMB.name} {V.shape}")


# =================================================================== umap
def cmd_umap(pca_dim=256, k=15):
    import faiss                      # never import torch in this process (OpenMP clash)
    import h5py
    import umap
    from sklearn.decomposition import PCA

    t0 = time.time()
    with h5py.File(H5, "r") as f:
        accessions = list(f.keys())
        X = np.vstack([f[a][()] for a in accessions]).astype(np.float32)
    print(f"loaded {X.shape[0]:,} x {X.shape[1]} in {time.time() - t0:.0f}s", flush=True)

    # The supplied h5 excludes all viral proteins (and ~11k others); embed_missing.py
    # computes those with the same verified pooling, so fold them in when present.
    if MISSING_EMB.exists() and MISSING_ACC.exists():
        extra_acc = MISSING_ACC.read_text().split()
        extra = np.load(MISSING_EMB).astype(np.float32)
        assert len(extra_acc) == len(extra), "missing-embedding accession/vector mismatch"
        new = [i for i, a in enumerate(extra_acc) if a not in set(accessions)]
        accessions += [extra_acc[i] for i in new]
        X = np.vstack([X, extra[new]])
        print(f"merged {len(new):,} separately embedded proteins -> {len(accessions):,} total",
              flush=True)

    q = np.load(ORF_EMB).astype(np.float32)
    names = ORF_IDS.read_text().split()
    X = np.vstack([X, q])

    # cosine geometry: L2-normalise, then euclidean ranks neighbours identically
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    t0 = time.time()
    X = np.ascontiguousarray(PCA(n_components=pca_dim, random_state=0).fit_transform(X))
    print(f"PCA -> {pca_dim}d in {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    index = faiss.IndexHNSWFlat(X.shape[1], 32)
    index.hnsw.efConstruction = 80
    index.add(X)
    index.hnsw.efSearch = 64
    dist, knn = index.search(X, k + 1)
    rows = np.arange(len(X))
    fix = knn[:, 0] != rows               # duplicates can outrank self
    knn[fix, 0], dist[fix, 0] = rows[fix], 0.0
    # Hand UMAP exactly n_neighbors columns: if it has to slice the array itself the
    # result is a non-contiguous view, which its numba kernels reject. Arrays must
    # also be C-contiguous for the same reason.
    knn = np.ascontiguousarray(knn[:, :k], dtype=np.int32)
    knn_dist = np.ascontiguousarray(np.sqrt(np.clip(dist[:, :k], 0, None)), dtype=np.float32)
    print(f"kNN graph in {time.time() - t0:.0f}s ({fix.mean():.1%} self-rank fixes)", flush=True)

    t0 = time.time()
    Y = umap.UMAP(n_neighbors=k, min_dist=0.1, metric="euclidean",
                  precomputed_knn=(knn, knn_dist, None), init="pca",
                  low_memory=True, n_jobs=-1, verbose=True).fit_transform(X)
    print(f"UMAP in {(time.time() - t0) / 60:.1f} min", flush=True)

    d = pd.DataFrame({"accession": accessions + names,
                      "is_query": [False] * len(accessions) + [True] * len(names),
                      "x": Y[:, 0], "y": Y[:, 1]})
    d = annotate(d)
    d.to_parquet(COORDS, index=False)
    print(f"Wrote {COORDS}")


def annotate(d):
    """Add taxonomic domain and InterPro homologous-superfamily labels."""
    import pyarrow.parquet as pq
    sys.path.insert(0, str(HERE))
    from annotate import domain

    meta = pq.read_table(PROJECT / "DS569k.parquet",
                         columns=["accession", "protein_name", "organism_name",
                                  "ncbi_taxonomy_phylum", "ncbi_taxonomy_class"]).to_pandas()
    meta = meta.set_index("accession")
    j = meta.reindex(d.accession)
    d["protein_name"] = j.protein_name.values
    d["organism_name"] = j.organism_name.values
    d["domain"] = [None if q else domain(p, c, o) for q, p, c, o in
                   zip(d.is_query, j.ncbi_taxonomy_phylum, j.ncbi_taxonomy_class,
                       j.organism_name)]

    sf_name = {}
    with open(ENTRY_LIST) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["ENTRY_TYPE"] == "Homologous_superfamily":
                sf_name[r["ENTRY_AC"]] = r["ENTRY_NAME"]
    want = set(d.accession[~d.is_query])
    sf = {}
    with gzip.open(XREFS, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Entry"] in want:
                ips = [x for x in r["InterPro"].split(";") if x in sf_name]
                if ips:
                    sf[r["Entry"]] = sf_name[ips[0]]
    d["superfamily"] = [None if q else sf.get(a) for q, a in zip(d.is_query, d.accession)]
    print(f"domain assigned: {d.domain.notna().mean():.1%}; "
          f"superfamily assigned: {d.superfamily.notna().sum():,} "
          f"({100 * d.superfamily.notna().sum() / (~d.is_query).sum():.1f}% of proteins)")
    return d


# =================================================================== plot
def cmd_plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    d = pd.read_parquet(COORDS)
    ref, q = d[~d.is_query], d[d.is_query]

    top = list(ref.superfamily.value_counts().head(N_SUPERFAMILIES).index)
    sf_colors = {s: matplotlib.colormaps["tab20"](i % 20) for i, s in enumerate(top)}
    OTHER = "#9a9994"
    # Panel D ranks superfamilies WITHIN the zoom window (zoom_top/zoom_colors, set below
    # once the window is known): the global top-20 covers only a few percent of that
    # region, which would otherwise leave the ORF neighbourhood almost entirely grey.

    # window around the ORFs for the zoom panels
    px = max(np.ptp(q.x) * 0.30, 1.5)
    py = max(np.ptp(q.y) * 0.30, 1.5)
    win = (q.x.min() - px, q.x.max() + px, q.y.min() - py, q.y.max() + py)
    inwin = ref.x.between(win[0], win[1]) & ref.y.between(win[2], win[3])

    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})
    fig, axes = plt.subplots(2, 2, figsize=(19, 17.5), facecolor=SURFACE)

    def base(ax, sub, title, zoom):
        ax.set_facecolor(SURFACE)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color(GRID); s.set_linewidth(0.8)
        ax.set_title(title, loc="left", fontsize=13, color=INK, pad=10)
        ax.scatter(q.x, q.y, s=190 if zoom else 110, c=ORF_COLOR, marker="*",
                   edgecolors="#ffffff", linewidths=1.2, zorder=6)
        if zoom:
            for _, r in q.iterrows():
                ax.annotate(r.accession, (r.x, r.y), xytext=(9, 9),
                            textcoords="offset points", fontsize=9.5, fontweight="bold",
                            color=INK, zorder=7,
                            bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none",
                                      alpha=0.85))
            ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3])

    def by_domain(ax, sub, title, zoom):
        base(ax, sub, title, zoom)
        s = 2.5 if zoom else 0.12
        miss = sub[sub.domain.isna()]
        if len(miss):
            ax.scatter(miss.x, miss.y, s=s, c=UNASSIGNED, linewidths=0, rasterized=True)
        for dom, col in DOMAIN_COLORS.items():
            m = sub[sub.domain == dom]
            ax.scatter(m.x, m.y, s=s, c=col, linewidths=0, alpha=0.65, rasterized=True)
        h = [Line2D([], [], marker="o", ls="", ms=8, mfc=c, mec="none",
                    label=f"{dom} ({(sub.domain == dom).sum():,})")
             for dom, c in DOMAIN_COLORS.items()]
        h.append(Line2D([], [], marker="*", ls="", ms=13, mfc=ORF_COLOR, mec="white",
                        label="query ORF"))
        ax.legend(handles=h, loc="upper right", frameon=True, facecolor=SURFACE,
                  edgecolor=GRID, fontsize=9.5, labelcolor=INK_2, framealpha=0.92,
                  title="Taxonomic domain", title_fontsize=10)

    def by_superfamily(ax, sub, title, zoom):
        base(ax, sub, title, zoom)
        s = 2.5 if zoom else 0.12
        palette, top_here = (zoom_colors, zoom_top) if zoom else (sf_colors, top)
        rest = sub[~sub.superfamily.isin(top_here)]
        na = rest[rest.superfamily.isna()]
        oth = rest[rest.superfamily.notna()]
        ax.scatter(na.x, na.y, s=s, c=UNASSIGNED, linewidths=0, rasterized=True)
        ax.scatter(oth.x, oth.y, s=s, c=OTHER, linewidths=0, alpha=0.5, rasterized=True)
        shown = [sf for sf in top_here if (sub.superfamily == sf).any()]
        for sf in shown:
            m = sub[sub.superfamily == sf]
            ax.scatter(m.x, m.y, s=s * 1.6, color=palette[sf], linewidths=0,
                       alpha=0.85, rasterized=True)
        h = [Line2D([], [], marker="o", ls="", ms=8, mfc=palette[sf], mec="none",
                    label=f"{sf[:44]}{'…' if len(sf) > 44 else ''} ({(sub.superfamily == sf).sum():,})")
             for sf in shown]
        h += [Line2D([], [], marker="o", ls="", ms=8, mfc=OTHER, mec="none",
                     label=f"other superfamily ({len(oth):,})"),
              Line2D([], [], marker="o", ls="", ms=8, mfc=UNASSIGNED, mec="none",
                     label=f"no superfamily assigned ({len(na):,})"),
              Line2D([], [], marker="*", ls="", ms=13, mfc=ORF_COLOR, mec="white",
                     label="query ORF")]
        ax.legend(handles=h, loc="upper right", frameon=True, facecolor=SURFACE,
                  edgecolor=GRID, fontsize=7.2, labelcolor=INK_2, framealpha=0.92,
                  title="Protein superfamily (InterPro homologous superfamily)",
                  title_fontsize=8.6, ncol=1)

    zoom_top = list(ref[inwin].superfamily.value_counts().head(N_ZOOM_SUPERFAMILIES).index)
    zoom_colors = {s: matplotlib.colormaps["tab20"](i % 20) for i, s in enumerate(zoom_top)}

    n = len(ref)
    by_domain(axes[0, 0], ref, f"A · All {n:,} Swiss-Prot proteins — taxonomic domain", False)
    by_superfamily(axes[0, 1], ref, f"B · All {n:,} Swiss-Prot proteins — protein superfamily", False)
    by_domain(axes[1, 0], ref[inwin], f"C · ORF neighbourhood ({inwin.sum():,} proteins) — taxonomic domain", True)
    by_superfamily(axes[1, 1], ref[inwin], f"D · ORF neighbourhood ({inwin.sum():,} proteins) — superfamily (ranked within this region)", True)

    for ax in (axes[0, 0], axes[0, 1]):
        ax.add_patch(plt.Rectangle((win[0], win[2]), win[1] - win[0], win[3] - win[2],
                                   fill=False, ec=INK_MUTED, lw=1.1, ls="--", zorder=5))

    fig.suptitle("UMAP of ESM-2 650M embeddings — Swiss-Prot and the query ORFs",
                 x=0.012, y=0.985, ha="left", fontsize=17, color=INK)
    fig.text(0.012, 0.963,
             f"{n:,} proteins (1280-d mean-pooled embeddings, PCA 256 -> FAISS kNN -> UMAP). "
             "Dashed box on A and B marks the region enlarged in C and D.",
             ha="left", fontsize=10.5, color=INK_2)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(FIG, dpi=200, facecolor=SURFACE)
    print(f"Wrote {FIG}")

    for _, r in q.iterrows():
        dist = np.hypot(ref.x - r.x, ref.y - r.y)
        nn = ref.loc[dist.nsmallest(20).index]
        sf = nn.superfamily.value_counts()
        dm = nn.domain.value_counts()
        print(f"{r.accession}: 20 nearest -> domains {dict(dm.head(3))}; "
              f"superfamilies {dict(sf.head(3)) if len(sf) else '(none assigned)'}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "verify":
        cmd_verify()
    elif cmd == "embed":
        cmd_embed()
    elif cmd == "umap":
        cmd_umap()
    elif cmd == "plot":
        cmd_plot()
    elif cmd == "all":
        cmd_embed(); cmd_umap(); cmd_plot()
    else:
        sys.exit(f"unknown command {cmd!r}")
