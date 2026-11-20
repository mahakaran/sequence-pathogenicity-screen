"""AlphaFold 3 structure prediction summary for the five ORFs.

Panel A  per-residue pLDDT along each ORF, all 5 seeds overlaid
Panel B  mean pLDDT vs MSA depth - the interpretive crux: the two confident
         models are the two with no homologs at all
Panel C  local helicity and chain extension, showing that the confident
         models are single helices rather than folds

Usage:  python af3_figure.py
"""

import glob
import re
import json
import math
import statistics as st
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
AF = PROJECT / "alphafold_orfs" / "extracted"
RES = HERE / "results"
FIG = RES / "af3_structures.png"

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
# sequential-by-magnitude for the pLDDT bands, categorical for the ORFs
ORF_C = {"ORF1": "#2a78d6", "ORF2": "#eda100", "ORF3": "#008300",
         "ORF4": "#e34948", "ORF5": "#7b4fa8"}
BANDS = [(0, 50, "#f2f1ee", "very low"), (50, 70, "#ffe9b8", "low"),
         (70, 90, "#c3ddf7", "confident"), (90, 100, "#9ec9f0", "very high")]

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

# ORF2 was submitted to AF3 as a two-chain job (chain A is ORF1's sequence),
# so ORF2's own prediction is chain B. Every other job is a single chain A.
QUERY_CHAIN = {"orf1": "A", "orf2": "B", "orf3": "A", "orf4": "A", "orf5": "A"}


def parse(path, chain):
    """Per-residue (aa, mean pLDDT) and CA coordinates for one chain."""
    res, ca = {}, {}
    for ln in open(path):
        if not ln.startswith("ATOM"):
            continue
        p = ln.split()
        if p[16] != chain:
            continue
        r = int(p[15])
        res.setdefault(r, [AA3.get(p[5], "X"), []])[1].append(float(p[14]))
        if p[3] == "CA":
            ca[r] = (float(p[10]), float(p[11]), float(p[12]))
    ks = sorted(res)
    return ks, "".join(res[k][0] for k in ks), [st.mean(res[k][1]) for k in ks], ca


def helix_frac(ca):
    """Fraction of i->i+4 CA distances in the alpha-helical range (5.0-7.5 A)."""
    ks = sorted(ca)
    ok = tot = 0
    for k in ks:
        if k + 4 in ca:
            tot += 1
            if 5.0 < math.dist(ca[k], ca[k + 4]) < 7.5:
                ok += 1
    return ok / tot if tot else 0.0


def msa_depth(orf, chain):
    """Aligned segments, including the query row (AF3's own 'depth')."""
    n = 0
    for f in glob.glob(f"{AF}/{orf}/msas/*chains_{chain.lower()}.a3m"):
        n = max(n, sum(1 for ln in open(f) if ln.startswith(">")))
    return n


def accession(header):
    """Bare accession from an a3m header, regardless of database spelling.

    AF3 writes the same protein two ways: 'UniRef90_G4YP12/40-107' in the
    unpaired file and 'tr|G4YP12|G4YP12_PHYSP/40-107' in the paired one. Without
    normalising, a union over both files counts every hit twice.
    """
    tok = header.lstrip(">").split()[0]
    tok = re.sub(r"/\d+-\d+.*$", "", tok)
    if "|" in tok:                       # tr|ACC|NAME  or  sp|ACC|NAME
        parts = tok.split("|")
        return parts[1] if len(parts) > 1 else tok
    return re.sub(r"^UniRef\d+_", "", tok)


def msa_proteins(orf, chain):
    """Distinct SOURCE PROTEINS behind those segments.

    AF3 reports a count of aligned regions, not of relatives: one repeat protein
    can contribute ten segments, and the paired file re-lists the unpaired hits.
    Collapsing to distinct accessions is what the depth actually means.
    """
    srcs = set()
    for f in glob.glob(f"{AF}/{orf}/msas/*chains_{chain.lower()}.a3m"):
        heads = [ln.strip() for ln in open(f) if ln.startswith(">")][1:]
        srcs |= {accession(h) for h in heads}
    return len(srcs)


def collect():
    out = {}
    for orf in sorted(QUERY_CHAIN):
        ch = QUERY_CHAIN[orf]
        seeds = []
        for m in sorted(glob.glob(f"{AF}/{orf}/*model_*.cif")):
            ks, seq, pl, ca = parse(m, ch)
            seeds.append((pl, ca))
        ks, seq, pl0, ca0 = parse(sorted(glob.glob(f"{AF}/{orf}/*model_0.cif"))[0], ch)
        s = json.load(open(glob.glob(f"{AF}/{orf}/*summary_confidences_0.json")[0]))
        pts = [ca0[k] for k in sorted(ca0)]
        cen = [st.mean(x) for x in zip(*pts)]
        out[orf.upper()] = {
            "seq": seq, "L": len(seq),
            "plddt": [s_[0] for s_ in seeds],
            "mean": st.mean(pl0),
            "ptm": s["ptm"], "disorder": s["fraction_disordered"],
            "msa": msa_depth(orf, ch),
            "msa_prot": msa_proteins(orf, ch),
            "helix": helix_frac(ca0),
            "rg": math.sqrt(st.mean([math.dist(p, cen) ** 2 for p in pts])),
            "e2e": math.dist(pts[0], pts[-1]),
            "n_templ": len(glob.glob(f"{AF}/{orf}/templates/*chains_{ch.lower()}.cif")),
        }
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    d = collect()
    names = list(d)

    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})
    fig = plt.figure(figsize=(16, 13.6), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.45, 1, 1],
                          top=0.878, bottom=0.075, left=0.068, right=0.982,
                          hspace=0.42, wspace=0.20)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])
    axD = fig.add_subplot(gs[2, :])

    def style(ax, grid="y"):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_2, labelsize=9)
        if grid:
            ax.grid(axis=grid, color=GRID, lw=0.7, alpha=0.8)
        ax.set_axisbelow(True)

    # ---- A: per-residue pLDDT, all seeds ---------------------------------
    style(axA)
    for lo, hi, c, _ in BANDS:
        axA.axhspan(lo, hi, color=c, zorder=0, alpha=0.75)
    for orf in names:
        for pl in d[orf]["plddt"]:
            axA.plot(range(1, len(pl) + 1), pl, color=ORF_C[orf], lw=1.0,
                     alpha=0.5, zorder=3)
        pl = d[orf]["plddt"][0]
        axA.plot(range(1, len(pl) + 1), pl, color=ORF_C[orf], lw=2.2,
                 alpha=0.95, zorder=4, label=f"{orf} ({d[orf]['L']} aa)")
        axA.annotate(orf, (len(pl), pl[-1]), xytext=(7, 0),
                     textcoords="offset points", fontsize=8.6, color=ORF_C[orf],
                     fontweight="bold", va="center", zorder=6)
    axA.axhline(70, color=INK_MUTED, lw=0.9, ls="--", zorder=5)
    axA.set_xlim(0, 103)
    axA.set_ylim(0, 100)
    axA.set_xlabel("residue position", fontsize=10.5, color=INK_2)
    axA.set_ylabel("pLDDT", fontsize=10.5, color=INK_2)
    axA.set_title("A · Per-residue confidence, all five seeds overlaid",
                  loc="left", fontsize=13, color=INK, pad=30)
    h = [Line2D([], [], color=ORF_C[o], lw=2.2, label=f"{o} ({d[o]['L']} aa)")
         for o in names]
    h += [Patch(fc=c, ec="none", label=lbl) for _, _, c, lbl in BANDS]
    axA.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 1.10),
               frameon=False, fontsize=9, labelcolor=INK_2, ncol=9,
               handletextpad=0.5, columnspacing=1.2)

    # ---- B: mean pLDDT vs MSA depth --------------------------------------
    style(axB)
    for orf in names:
        x = max(d[orf]["msa_prot"], 1)
        axB.scatter([x], [d[orf]["mean"]], s=230, color=ORF_C[orf],
                    edgecolors="white", linewidths=1.6, zorder=5)
        # ORF3 and ORF5 both sit at MSA depth 1 and within 6 pLDDT of each other
        OFFB = {"ORF5": (12, 4), "ORF3": (12, -20), "ORF1": (12, -20),
                "ORF4": (-12, 6), "ORF2": (10, 6)}
        off = OFFB.get(orf, (10, 6))
        axB.annotate(f"{orf}\n{d[orf]['mean']:.0f}", (x, d[orf]["mean"]),
                     xytext=off, textcoords="offset points", fontsize=8.4,
                     color=INK, fontweight="bold", zorder=6,
                     bbox=dict(boxstyle="round,pad=0.16", fc=SURFACE, ec="none",
                               alpha=0.88))
    axB.axhline(70, color=INK_MUTED, lw=0.9, ls="--", zorder=2)
    axB.text(40, 71.5, "confident", fontsize=8, color=INK_MUTED, va="bottom",
             ha="right")
    axB.set_xscale("log")
    axB.set_xlim(0.8, 60)
    axB.set_ylim(0, 105)
    axB.set_xticks([1, 2, 5, 10, 20, 50])
    axB.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axB.set_xlabel("distinct source proteins in the MSA (1 = query only)",
                   fontsize=10.5, color=INK_2)
    axB.set_ylabel("mean pLDDT", fontsize=10.5, color=INK_2)
    axB.set_title("B · Confidence against homolog count", loc="left",
                  fontsize=13, color=INK, pad=10)

    # ---- C: helicity vs extension ----------------------------------------
    style(axC)
    for orf in names:
        x = d[orf]["helix"]
        y = d[orf]["e2e"] / d[orf]["L"]
        axC.scatter([x], [y], s=230, color=ORF_C[orf], edgecolors="white",
                    linewidths=1.6, zorder=5)
        axC.annotate(orf, (x, y), xytext=(10, 5), textcoords="offset points",
                     fontsize=8.4, color=INK, fontweight="bold", zorder=6,
                     bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none",
                               alpha=0.88))
    axC.axhline(1.35, color=INK_MUTED, lw=0.9, ls="--", zorder=2)
    axC.text(0.02, 1.38, "extended rod (≈1.5 Å per residue)", fontsize=8,
             color=INK_MUTED, va="bottom")
    axC.set_xlim(-0.05, 1.12)
    axC.set_ylim(0, 1.75)
    axC.set_xlabel("fraction of i→i+4 Cα distances in the α-helical range",
                   fontsize=10.5, color=INK_2)
    axC.set_ylabel("end-to-end distance / residue (Å)", fontsize=10.5, color=INK_2)
    axC.set_title("C · Helix content against chain extension", loc="left",
                  fontsize=13, color=INK, pad=10)

    # ---- D: metric table --------------------------------------------------
    axD.axis("off")
    axD.set_facecolor(SURFACE)
    # pin the table's own coordinate system: without explicit limits the text
    # below lands on whatever autoscale picks and overruns the panels above
    axD.set_xlim(0, 1)
    axD.set_ylim(0, 1)
    cols = ["ORF", "Len", "Proteins (segs)", "Mean pLDDT", "Residues >70", "pTM",
            "Fraction disordered", "Helical frac", "Interpretation"]
    wid = [0.055, 0.042, 0.082, 0.085, 0.092, 0.050, 0.120, 0.082, 0.392]
    xs, acc = [], 0.0
    for w in wid:
        xs.append(acc)
        acc += w
    yy = 0.86
    for x, w, c in zip(xs, wid, cols):
        axD.text(x, yy, c, fontsize=9, color=INK, fontweight="bold", va="top")
    axD.plot([0, 1], [yy - 0.10, yy - 0.10], color=INK_2, lw=1.0)
    yy -= 0.20
    INTERP = {
        "ORF1": "no fold; low confidence throughout",
        "ORF2": "no fold; lowest confidence of the five",
        "ORF3": "single short helix, no tertiary structure",
        "ORF4": "no fold; low confidence throughout",
        "ORF5": "single long hydrophobic helix (TM-like)",
    }
    for orf in names:
        v = d[orf]
        n70 = sum(1 for p in v["plddt"][0] if p > 70)
        cells = [orf, str(v["L"]), f"{v['msa_prot']} ({v['msa']})", f"{v['mean']:.1f}",
                 f"{n70}/{v['L']}", f"{v['ptm']:.2f}", f"{v['disorder']:.2f}",
                 f"{v['helix']:.2f}", INTERP[orf]]
        for i, (x, cell) in enumerate(zip(xs, cells)):
            axD.text(x, yy, cell, fontsize=8.8,
                     color=ORF_C[orf] if i == 0 else INK_2,
                     fontweight="bold" if i == 0 else "normal", va="top")
        yy -= 0.155
    axD.set_title("D · AlphaFold 3 confidence metrics and structural character",
                  loc="left", fontsize=13, color=INK, pad=10)

    fig.suptitle("AlphaFold 3 predictions for the five ORFs",
                 x=0.068, y=0.978, ha="left", fontsize=17, color=INK)
    fig.text(0.068, 0.950,
             "AlphaFold 3 server, 5 seeds per ORF. Confidence is reproducible across seeds "
             "(mean-pLDDT spread ≤ 2.7 for every ORF).\nORF2 was submitted as a two-chain job "
             "in which chain A repeats ORF1's sequence; the values here are ORF2's own chain B.",
             ha="left", fontsize=10.3, color=INK_2, va="top")
    fig.text(0.068, 0.030,
             "The two models that reach confident pLDDT are the two for which "
             "AlphaFold found no homologs at all (MSA depth 1 = the query alone), and "
             "both are single helices with no tertiary\nstructure — AlphaFold's own "
             "fraction_disordered is 1.00 for each. High local confidence in a lone "
             "helix is not evidence of a fold, and carries no functional information.",
             fontsize=8.6, color=INK_MUTED, ha="left")

    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=200, facecolor=SURFACE)
    print(f"Wrote {FIG}")
    for orf in names:
        v = d[orf]
        print(f"  {orf}: L={v['L']} msa={v['msa']} templ={v['n_templ']} "
              f"mean_pLDDT={v['mean']:.1f} pTM={v['ptm']:.2f} "
              f"disorder={v['disorder']:.2f} helix={v['helix']:.2f} "
              f"Rg={v['rg']:.1f} e2e={v['e2e']:.1f}")


if __name__ == "__main__":
    main()
