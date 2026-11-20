"""Render the ORF-vs-flagged-neighbour protein alignments as viewable figures.

One figure per ORF, showing every flagged neighbour's local Smith-Waterman
alignment (BLOSUM62, gaps -11/-1) in full, ordered by empirical p-value.

Colouring marks what the eye should judge:
  identical residue   bold ink on a tinted cell
  conserved (BLOSUM>0) mid ink
  everything else     muted

A run of >=4 consecutive identities is boxed: that is the shape real homology
makes (a conserved motif), and its absence is the finding here.

Usage:  python alignment_figures.py [--max-per-orf 12]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
RES = HERE / "results"
FIGDIR = RES / "alignments"

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
IDENT_BG = "#cfe0f5"       # tint behind an identical residue
MOTIF = "#2a78d6"          # box around a >=4-run of identities
ORF_COLOR = "#eb6834"
WRAP = 60                  # residues per line
BL = substitution_matrices.load("BLOSUM62")


def conserved(a, b):
    if a in "- " or b in "- ":
        return False
    try:
        return BL[a, b] > 0
    except (KeyError, IndexError):
        return False


def runs_of_identity(q, s, min_len=4):
    """Spans (start, end) where query and subject are identical for >=min_len."""
    out, start = [], None
    for i, (a, b) in enumerate(zip(q, s)):
        if a == b and a != "-":
            start = i if start is None else start
        else:
            if start is not None and i - start >= min_len:
                out.append((start, i))
            start = None
    if start is not None and len(q) - start >= min_len:
        out.append((start, len(q)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-per-orf", type=int, default=12,
                    help="show this many alignments per ORF, lowest p-value first")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    d = pd.read_parquet(RES / "flagged_alignments.parquet")
    meta = pq.read_table(PROJECT / "DS569k.parquet",
                         columns=["accession", "sequence"]).to_pandas().set_index("accession")
    seqs = {}
    for blk in (PROJECT / "screening-exercise-orfs.fasta").read_text().split(">"):
        if blk.strip():
            h, *rest = blk.strip().split("\n")
            seqs[h.split("|")[-1].split(":")[0]] = "".join(rest).rstrip("*").upper()

    al = Align.PairwiseAligner(scoring="blastp")
    al.mode = "local"
    FIGDIR.mkdir(parents=True, exist_ok=True)
    bonferroni = 0.05 / len(d)

    for orf, grp in d.groupby("orf"):
        grp = grp.sort_values("p_empirical").head(args.max_per_orf)
        q_full = "".join(c for c in seqs[orf] if c in "ACDEFGHIKLMNPQRSTVWY")

        blocks = []
        for _, r in grp.iterrows():
            s = "".join(c for c in str(meta.loc[r.accession, "sequence"]).upper()
                        if c in "ACDEFGHIKLMNPQRSTVWY")
            aln = al.align(q_full, s)[0]
            # Take the aligned strings from the alignment object itself; parsing the
            # printed form drags in Biopython's coordinate columns and padding spaces.
            a_q, a_s = aln[0], aln[1]
            n = min(len(a_q), len(a_s))
            blocks.append((r, a_q[:n], a_s[:n]))

        rows_needed = sum(int(np.ceil(len(b[1]) / WRAP)) for b in blocks)
        fig_h = 1.5 + 1.25 * len(blocks) + 0.46 * rows_needed
        fig, ax = plt.subplots(figsize=(15.5, fig_h), facecolor=SURFACE)
        ax.set_facecolor(SURFACE)
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, fig_h)
        y = fig_h - 0.35

        ax.text(0, y, f"{orf} — protein–protein alignments vs flagged neighbours",
                fontsize=15, color=INK, fontweight="bold", va="top")
        y -= 0.34
        ax.text(0, y,
                f"{len(q_full)} aa · local Smith-Waterman, BLOSUM62, gaps −11/−1 · "
                f"ordered by empirical p (2,000 composition-preserving shuffles) · "
                f"Bonferroni threshold for 92 tests: p < {bonferroni:.5f}",
                fontsize=9, color=INK_2, va="top")
        y -= 0.30
        ax.text(0, y, "ORF sequence:  " + q_full, fontsize=7.6, color=INK_MUTED,
                family="monospace", va="top")
        y -= 0.52

        cw = 0.0118        # width of one monospace column in axes units
        for r, a_q, a_s in blocks:
            sig = "n.s." if r.p_empirical >= bonferroni else "SIGNIFICANT"
            flags = []
            if r.pathogen == "yes":
                flags.append("pathogen source")
            if r.pathogen == "venom":
                flags.append("venomous source")
            if r.vf:
                flags.append(f"VFDB: {r.vf}")
            ax.text(0, y, f"#{int(r['rank'])}  {r.accession}  {str(r.protein_name)[:52]}",
                    fontsize=9.6, color=INK, fontweight="bold", va="top")
            y -= 0.235
            ax.text(0, y, f"{str(r.organism)[:46]}   ·   " + " · ".join(flags),
                    fontsize=8.1, color=INK_2, va="top")
            y -= 0.225
            ax.text(0, y,
                    f"SW {r.sw_score:.0f}   z {r.z:+.2f}   p {r.p_empirical:.4f} ({sig})   "
                    f"identity {r.pct_identity:.0f}% over {int(r.aln_length)} aligned residues "
                    f"({100 * r.aln_length / len(q_full):.0f}% of {orf})   "
                    f"cosine dist {r.cosine_distance:.4f}",
                    fontsize=8.1, color=INK_2, va="top")
            y -= 0.34

            for off in range(0, len(a_q), WRAP):
                cq, cs = a_q[off:off + WRAP], a_s[off:off + WRAP]
                spans = runs_of_identity(cq, cs)
                for st, en in spans:                       # motif box behind the run
                    ax.add_patch(Rectangle((0.058 + st * cw, y - 0.40),
                                           (en - st) * cw, 0.42,
                                           fill=False, ec=MOTIF, lw=1.2, zorder=3))
                ax.text(0.0, y, f"{orf}", fontsize=7.4, color=INK_MUTED,
                        family="monospace", va="top")
                for i, (a, b) in enumerate(zip(cq, cs)):
                    x = 0.058 + i * cw
                    if a == b and a != "-":
                        ax.add_patch(Rectangle((x - 0.0012, y - 0.135), cw, 0.155,
                                               color=IDENT_BG, lw=0, zorder=1))
                    ax.text(x, y, a, fontsize=7.4, family="monospace", va="top", zorder=2,
                            color=INK if a == b and a != "-" else
                            (INK_2 if conserved(a, b) else INK_MUTED),
                            fontweight="bold" if a == b and a != "-" else "normal")
                y -= 0.175
                ax.text(0.058, y, "".join("|" if a == b and a != "-" else
                                          (":" if conserved(a, b) else " ")
                                          for a, b in zip(cq, cs)),
                        fontsize=7.4, family="monospace", color=INK_MUTED, va="top")
                y -= 0.175
                ax.text(0.0, y, "neigh", fontsize=7.4, color=INK_MUTED,
                        family="monospace", va="top")
                for i, (a, b) in enumerate(zip(cq, cs)):
                    ax.text(0.058 + i * cw, y, b, fontsize=7.4, family="monospace",
                            va="top", zorder=2,
                            color=INK if a == b and a != "-" else
                            (INK_2 if conserved(a, b) else INK_MUTED),
                            fontweight="bold" if a == b and a != "-" else "normal")
                y -= 0.30
            y -= 0.18
            ax.plot([0, 1], [y, y], color=GRID, lw=0.7)
            y -= 0.30

        ax.text(0, y, "Boxed = run of ≥4 consecutive identical residues, the shape a conserved "
                      "motif makes. Bold+tinted = identical; mid-tone = conserved (BLOSUM62 > 0).",
                fontsize=8, color=INK_MUTED, va="top")
        out = FIGDIR / f"{orf}_alignments.png"
        fig.savefig(out, dpi=185, facecolor=SURFACE, bbox_inches="tight")
        plt.close(fig)
        n_motif = sum(len(runs_of_identity(b[1], b[2])) for b in blocks)
        print(f"  {out.name}: {len(blocks)} alignments, {n_motif} runs of ≥4 identities")


if __name__ == "__main__":
    main()
