"""Infographic cover page: the question, the evidence, and the finding at a glance.

Drawn at A4 proportions so it can be dropped in as page 1 of the report. Every
number here is taken from the report body; nothing is restated from memory.

Layout: a top-down cursor in figure-fraction units. Each block is placed at the
cursor and then advances it by its own measured height plus an explicit gap, so
blocks cannot collide no matter how their contents change. Vertical positions are
never hardcoded - an earlier version of this script tuned ~20 individual
y-constants by hand, and every text change reintroduced an overlap elsewhere.

Usage:  python cover_figure.py
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "cover.png"

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
NEG = "#2a78d6"      # a negative / no-hit result
POS = "#eda100"      # the one positive finding
WARN = "#e34948"
GOOD = "#008300"
AMBER_INK, RED_INK = "#8a6400", "#9c2f2f"

W, H = 8.27, 11.69   # A4 in inches
L, R = 0.072, 0.928  # left / right text margins, in figure fraction


def lh(fontsize, leading=1.0):
    """Height of one text line, as a fraction of figure height."""
    return fontsize * leading / 72.0 / H


QUERY = ("ATGAATCAACATAATACCCAAATCAATAAATTTATCTTTTTTAGTAGCTTCAGAATCATCAGTACCACCCCAATTCAATTT"
         "ATCAACAATAATAGGACCAGCAACAGCAATAGCAACATATTCATCATTAGAAGCAGGAGAATTCAAAACAATACCAACAGC"
         "TTTTTTATTATCAGCAGTATTCCAAGTTTTACAAACAGCACCATTAGAATTACCTTTTTCACAATCACCAGCATGCAAATT"
         "ATCAATACCAGAATTTTTATGAGATCTAATATGAACCAAACCCAT")

EVIDENCE = [
    ("BLASTn", "NCBI nt, E ≤ 1000", "0 hits"),
    ("DIAMOND blastp", "Swiss-Prot, 569,192 proteins", "0 hits"),
    ("jackhmmer ×5", "Swiss-Prot, profile search", "0 recruited"),
    ("AlphaFold 3 MSA", "AF3 databases", "none survive"),
    ("Kraken2 ×2 DBs", "Standard-8 + PlusPF-8", "0 / 254 k-mers"),
    ("DIAMOND vs VFDB", "4,755 virulence factors", "0 hits"),
]

ORFS = [
    ("ORF1", 96, 41.6, "no fold"),
    ("ORF2", 86, 28.7, "no fold"),
    ("ORF3", 25, 83.6, "lone helix"),
    ("ORF4", 96, 41.2, "no fold"),
    ("ORF5", 47, 89.5, "lone helix"),
]

PLL_MARKS = [(-3.123, "uniform random", INK_MUTED),
             (-2.895, "shuffled", INK_MUTED),
             (-2.688, "ORF residues shuffled", INK_MUTED),
             (-1.730, "real Swiss-Prot proteins", GOOD)]


def track(t, n=1):
    """Letter-spaced small caps: matplotlib Text has no letterspacing property."""
    return (" " * n).join(t)


def wrap(fig, s, fontsize, max_w, weight="normal"):
    """Greedy word-wrap to a width in figure fractions, measuring the real
    rendered text. A hand-placed newline does not know where a panel's right
    edge is, so a long string silently ran past it."""
    r = fig.canvas.get_renderer()

    def width(t):
        probe = fig.text(0, 0, t, fontsize=fontsize, fontweight=weight)
        w = probe.get_window_extent(renderer=r).width / fig.bbox.width
        probe.remove()
        return w

    out, line = [], ""
    for word in s.split():
        trial = word if not line else line + " " + word
        if line and width(trial) > max_w:
            out.append(line)
            line = word
        else:
            line = trial
    if line:
        out.append(line)
    return "\n".join(out)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Rectangle

    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_facecolor(SURFACE)

    y = 0.968                       # the cursor: top edge of the next block

    def text(s, size, color=INK, x=L, weight="normal", style="normal",
             ha="left", leading=1.30, zorder=5):
        """Draw a (possibly multi-line) string with its TOP at the cursor."""
        nonlocal y
        ax.text(x, y, s, fontsize=size, color=color, va="top", ha=ha,
                fontweight=weight, style=style, linespacing=leading, zorder=zorder)
        y -= lh(size, leading) * s.count("\n") + lh(size)

    def gap(frac):
        nonlocal y
        y -= frac

    def band(height, fc, accent=None, pad=0.020):
        """Rounded panel of `height` hanging from the cursor. Returns its top edge
        and inner text x, and leaves the cursor at the panel's bottom."""
        nonlocal y
        top = y
        ax.add_patch(FancyBboxPatch((L, top - height), R - L, height,
                                    boxstyle="round,pad=0,rounding_size=0.008",
                                    fc=fc, ec="none", zorder=1))
        if accent:
            ax.add_patch(Rectangle((L, top - height), 0.0045, height, fc=accent,
                                   ec="none", zorder=3))
        y = top - height
        return top, L + pad

    # ---------------- header ------------------------------------------------
    text(track("SEQUENCE SCREENING REPORT"), 8.6, INK_MUTED, weight="bold")
    gap(0.010)
    text("Is this sequence\npathogenic?", 26, INK, weight="bold", leading=1.14)
    gap(0.008)
    text("An unknown 288 bp nucleotide sequence · 5 open reading frames · "
         "no supplied provenance", 9.6, INK_2)
    gap(0.010)
    ax.plot([L, R], [y, y], color=GRID, lw=1.1, zorder=4)
    gap(0.010)

    # ---------------- the query --------------------------------------------
    text(track("THE QUERY"), 7.2, INK_MUTED, weight="bold")
    gap(0.005)
    for i in range(2):
        ax.text(L, y, QUERY[i * 144:(i + 1) * 144], fontsize=5.6, color=INK_MUTED,
                family="monospace", va="top", zorder=5)
        y -= lh(5.6, 1.60)
    gap(0.014)

    # ---------------- the answer -------------------------------------------
    top, tx = band(0.094, "#eef3fb")
    yy = top - 0.018
    ax.text(tx, yy, track("THE FINDING"), fontsize=7.4, color=NEG, va="top",
            fontweight="bold", zorder=5)
    yy -= lh(7.4) + 0.009
    ax.text(tx, yy, "No evidence of pathogenicity.\nNo evidence of natural origin either.",
            fontsize=14.0, color=INK, va="top", fontweight="bold", linespacing=1.24,
            zorder=5)
    yy -= lh(14.0, 1.24) + lh(14.0) + 0.009
    ax.text(tx, yy, "Most probable interpretation: a designed sequence, back-translated "
                    "from protein, with no natural homolog.",
            fontsize=9.0, color=INK_2, va="top", zorder=5)
    gap(0.015)

    # ---------------- evidence grid ----------------------------------------
    ax.text(R, y, "every database searched is homology-based", fontsize=7.6,
            color=INK_MUTED, va="top", ha="right", style="italic", zorder=5)
    text(track("SIX SEARCHES, NO HOMOLOG"), 7.4, INK_MUTED, weight="bold")
    gap(0.013)

    colw = (R - L) / 3
    rowgap = 0.0105
    # Card height is derived from the three stacked lines rather than guessed: a
    # hardcoded height shorter than its contents pushed the result line through
    # the bottom border.
    CARD_PAD_T, CARD_PAD_B = 0.0088, 0.0088
    CARD_LINES = ((8.8, 0.0040), (6.9, 0.0058), (9.2, 0.0))   # (fontsize, gap below)
    cardh = CARD_PAD_T + sum(lh(f) + g for f, g in CARD_LINES) + CARD_PAD_B
    for i, (name, db, res) in enumerate(EVIDENCE):
        cx = L + (i % 3) * colw
        ctop = y - (i // 3) * (cardh + rowgap)
        ax.add_patch(FancyBboxPatch((cx, ctop - cardh), colw - 0.014, cardh,
                                    boxstyle="round,pad=0,rounding_size=0.006",
                                    fc=SURFACE, ec=GRID, lw=0.9, zorder=2))
        ax.add_patch(Rectangle((cx, ctop - cardh), 0.0045, cardh, fc=NEG,
                               ec="none", zorder=3))
        cy = ctop - CARD_PAD_T
        for (s_, (fs, gapb), col, wt) in zip(
                (name, db, res), CARD_LINES,
                (INK, INK_MUTED, NEG), ("bold", "normal", "bold")):
            ax.text(cx + 0.016, cy, s_, fontsize=fs, color=col, va="top",
                    fontweight=wt, zorder=4)
            cy -= lh(fs) + gapb
    y -= 2 * cardh + rowgap
    gap(0.009)

    # ---------------- the one positive finding ------------------------------
    codon_body = wrap(fig,
                      "Back-translating ORF4's peptide with that 19-entry table reproduces "
                      "the 288 bp exactly. The three codons whose reverse complements are "
                      "stop codons (CTA, TCA, TTA) are absent entirely \u2014 which is what "
                      "keeps both strands open.", 8.3, (R - L) - 2 * 0.020)
    codon_h = (0.016 + lh(7.4) + 0.011 + lh(10.4) + 0.010
               + lh(8.3, 1.45) * codon_body.count("\n") + lh(8.3) + 0.013)
    top, tx = band(codon_h, "#fdf6e6", accent=POS)
    yy = top - 0.016
    ax.text(tx, yy, track("THE ONE POSITIVE FINDING — CODON USAGE"),
            fontsize=7.4, color=AMBER_INK, va="top", fontweight="bold", zorder=5)
    ax.text(R - 0.020, yy, "no natural gene behaves this way", fontsize=8.4,
            color=AMBER_INK, va="top", ha="right", fontweight="bold", zorder=5)
    yy -= lh(7.4) + 0.011
    ax.text(tx, yy, "ENC = 20.0, the theoretical floor: exactly one codon per amino acid.",
            fontsize=10.4, color=INK, va="top", fontweight="bold", zorder=5)
    yy -= lh(10.4) + 0.010
    ax.text(tx, yy,
            codon_body,
            fontsize=8.3, color=INK_2, va="top", linespacing=1.45, zorder=5)
    gap(0.012)

    # ---------------- structure strip ---------------------------------------
    text(track("WHAT THE FIVE PEPTIDES LOOK LIKE"), 7.4, INK_MUTED, weight="bold")
    gap(0.004)
    text("AlphaFold 3 mean pLDDT (0–100); dashed line = 70, the confidence threshold",
         7.2, INK_MUTED)
    gap(0.015)

    bar_top = y
    h_full = 0.062
    bw = (R - L - 0.30) / 5
    for i, (name, L_aa, pl, verdict) in enumerate(ORFS):
        bx = L + i * bw
        c = GOOD if pl > 70 else NEG
        base = bar_top - h_full
        ax.add_patch(Rectangle((bx, base), bw - 0.014, h_full, fc="#efeeea",
                               ec="none", zorder=2))
        ax.add_patch(Rectangle((bx, base), bw - 0.014, h_full * pl / 100, fc=c,
                               ec="none", alpha=0.85, zorder=3))
        ax.plot([bx, bx + bw - 0.014], [base + h_full * 0.70] * 2,
                color=INK_MUTED, lw=0.7, ls=(0, (2, 2)), zorder=4)
        # value above the fill, never inside it: at pLDDT 29 the fill is too
        # shallow to carry readable text
        ax.text(bx + (bw - 0.014) / 2, base + h_full * pl / 100 + 0.0035,
                f"{pl:.0f}", fontsize=8.8, color=c, ha="center", va="bottom",
                fontweight="bold", zorder=5)
        ly = base - 0.010
        for s, fs, col, wt in ((name, 8.4, INK, "bold"),
                               (f"{L_aa} aa", 7.0, INK_MUTED, "normal"),
                               (verdict, 6.9, c, "bold")):
            ax.text(bx + (bw - 0.014) / 2, ly, s, fontsize=fs, color=col,
                    ha="center", va="top", fontweight=wt, zorder=5)
            ly -= lh(fs) + 0.005

    strip_h = h_full + 0.010 + lh(8.4) + lh(7.0) + lh(6.9) + 3 * 0.005 + 0.004
    cx = L + 5 * bw + 0.012
    ax.add_patch(FancyBboxPatch((cx, bar_top - strip_h), R - cx, strip_h,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                fc=SURFACE, ec=GRID, lw=0.9, zorder=2))
    cyy = bar_top - 0.014
    ax.text(cx + 0.015, cyy, "The two confident models", fontsize=8.4, color=INK,
            va="top", fontweight="bold", zorder=3)
    cyy -= lh(8.4) + 0.009
    ax.text(cx + 0.015, cyy,
            "ORF3 and ORF5 are the only ones\nAlphaFold models confidently — and\nthe only "
            "two with an MSA of depth 1,\ni.e. no homolog found at all.",
            fontsize=7.3, color=INK_2, va="top", linespacing=1.5, zorder=3)
    cyy -= lh(7.3, 1.5) * 3 + lh(7.3) + 0.010
    ax.text(cx + 0.015, cyy,
            "Both are lone helices, not folds:\nfraction_disordered = 1.00.",
            fontsize=7.3, color=INK_2, va="top", linespacing=1.5, zorder=3)
    y = bar_top - strip_h
    gap(0.012)

    # ---------------- protein-likeness --------------------------------------
    text(track("ARE THEY PROTEIN-LIKE?  ESM-2 PSEUDO-LIKELIHOOD"), 7.4,
         INK_MUTED, weight="bold")
    gap(0.015)

    sx0, sx1 = L, R
    lo, hi = -3.26, -1.58

    def sx(v):
        return sx0 + (v - lo) / (hi - lo) * (sx1 - sx0)

    ax.text(sx(-2.565), y, "the five ORFs  −2.39 to −2.74", fontsize=8.0,
            color=WARN, ha="center", va="top", fontweight="bold", zorder=5)
    y -= lh(8.0) + 0.008

    track_h = 0.012
    ax.add_patch(Rectangle((sx0, y - track_h), sx1 - sx0, track_h, fc="#efeeea",
                           ec="none", zorder=2))
    ax.add_patch(Rectangle((sx(-2.74), y - track_h - 0.004),
                           sx(-2.39) - sx(-2.74), track_h + 0.008,
                           fc=WARN, ec="none", alpha=0.85, zorder=3))
    for v, lbl, c in PLL_MARKS:
        ax.plot([sx(v), sx(v)], [y - track_h - 0.005, y + 0.002], color=c,
                lw=1.6, zorder=4)
    y -= track_h + 0.010
    # first and last tick labels are clamped to the margins so they cannot run
    # off the page edge
    for v, lbl, c in PLL_MARKS:
        if v == PLL_MARKS[0][0]:
            ha, xpos = "left", sx0
        elif v == PLL_MARKS[-1][0]:
            ha, xpos = "right", sx1
        else:
            ha, xpos = "center", sx(v)
        ax.text(xpos, y, lbl, fontsize=6.8, color=c, ha=ha, va="top", zorder=5)
    y -= lh(6.8) + 0.010
    text(wrap(fig,
              "The ORFs sit inside the negative-control band. Against their own "
              "composition-matched shuffles ORF1 and ORF5 sit at the 50th percentile — "
              "their own residues in random order are exactly as protein-like as the real "
              "sequence.", 7.6, R - L),
         7.6, INK_2, leading=1.5)
    gap(0.011)

    # ---------------- the biosecurity point ---------------------------------
    red_body = wrap(fig,
                    "Every resource used here detects hazard by resemblance to catalogued "
                    "examples. A designed sequence with no natural homolog is invisible to "
                    "all of them \u2014 whether it is benign or not. The absence of hits is "
                    "uninformative about function.", 8.5, (R - L) - 2 * 0.020)
    red_h = 0.014 + lh(7.4) + 0.010 + lh(8.5, 1.5) * red_body.count("\n") + lh(8.5) + 0.013
    top, tx = band(red_h, "#fbeeee", accent=WARN)
    yy = top - 0.014
    ax.text(tx, yy, track("WHY A CLEAN SCREEN IS NOT A CLEARANCE"), fontsize=7.4,
            color=RED_INK, va="top", fontweight="bold", zorder=5)
    yy -= lh(7.4) + 0.010
    ax.text(tx, yy,
            red_body,
            fontsize=8.5, color=INK, va="top", linespacing=1.5, zorder=5)
    gap(0.014)

    # ---------------- footer ------------------------------------------------
    ax.plot([L, R], [y, y], color=GRID, lw=1.0, zorder=4)
    gap(0.011)
    text(wrap(fig,
              "Full methods, figures, limitations and conclusion follow.  Analyses: "
              "homology search · taxonomic classification · codon usage · protein "
              "language-model embedding · pairwise alignment · structure prediction "
              "· pseudo-likelihood modelling.", 7.2, R - L),
         7.2, INK_MUTED, leading=1.5)

    # The cursor must finish above the page edge. Anything below ~0.008 means a
    # block has been added without reclaiming space elsewhere, and the layout
    # will start colliding at the bottom.
    if y < 0.004:
        print(f"WARNING: content overruns the page bottom (cursor at {y:+.4f})")
    else:
        print(f"bottom margin clear: cursor ended at {y:+.4f}")

    fig.savefig(OUT, dpi=240, facecolor=SURFACE)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
