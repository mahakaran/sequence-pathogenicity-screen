"""Methods-overview figure: what each tool asks, and at what level of abstraction.

The argument the figure makes is that the eight analyses are not eight ways of
asking the same question. They sit on a ladder of abstraction - exact nucleotide
matching at the bottom, learned representations of protein space at the top - and
a sequence can be invisible to everything at one level while still being
detectable at another. The right-hand column states what a positive result would
have looked like, which is what makes the uniform negative interpretable.

Layout uses the same measured-cursor discipline as cover_figure.py: nothing is
positioned by a hardcoded y, so blocks cannot collide.

Usage:  python approach_figure.py
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "approach.png"

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
NUC, PROT, MODEL, STRUCT, ASSESS = "#2a78d6", "#008300", "#7b4fa8", "#eda100", "#e34948"

W, H = 12.4, 5.60
L, R = 0.050, 0.950


def lh(fontsize, leading=1.0):
    return fontsize * leading / 72.0 / H


# (level label, colour, [(tool, what it does, what a hit would have meant)])
LEVELS = [
    ("NUCLEOTIDE\nthe letters\nthemselves", NUC, [
        ("BLASTn",
         "Seeds on exact short word matches (here, 7 nt), then extends them into "
         "gapped local alignments and scores each against what chance would give.",
         "a related gene anywhere\nin GenBank"),
        ("Kraken2",
         "Breaks the sequence into every overlapping 35-mer and looks each one up "
         "in a hash of genomes, assigning the lowest common ancestor of the matches.",
         "the organism, or at least\nthe clade, it came from"),
        ("Composition & codon usage",
         "GC content, CpG frequency and which synonymous codon is chosen for each "
         "amino acid — signatures a genome imprints on its genes over evolutionary time.",
         "a host genome whose\nsignature it shares"),
    ]),
    ("PROTEIN\nthe translated\npeptides", PROT, [
        ("DIAMOND blastp",
         "Same logic as BLAST but on amino acids, scoring substitutions by chemical "
         "similarity (BLOSUM62), so divergent-but-related proteins still align.",
         "a homologous protein\nin Swiss-Prot or VFDB"),
        ("jackhmmer",
         "Iterative profile search: builds a position-specific model from whatever it "
         "finds, re-searches with it, and repeats — detecting homology too remote to pair.",
         "a protein family, found\nthrough intermediates"),
        ("Pairwise alignment vs flagged hits",
         "Full Smith-Waterman against every neighbour from a pathogen or virulence "
         "database, calibrated against shuffles of the query that preserve composition.",
         "real similarity to a toxin,\nnot a compositional echo"),
    ]),
    ("LEARNED &\nSTRUCTURAL\nbeyond detectable\nhomology", MODEL, [
        ("ESM-2 embedding",
         "A language model trained on 65M proteins maps each sequence to a vector; "
         "proteins with similar structure and function land near one another.",
         "a functional neighbourhood\nwhere sequence search failed"),
        ("ESM-2 pseudo-likelihood",
         "Masks each residue and asks how well the model predicts it — scoring how "
         "protein-like the sequence is, independent of any specific relative.",
         "a sequence evolution\nplausibly produced"),
        ("AlphaFold 3",
         "Predicts 3D structure from sequence and its own alignment search. Structure "
         "outlives sequence similarity, so a fold can betray an origin the letters hide.",
         "a fold matching a known\ntoxin or virulence scaffold"),
    ]),
]


def track(t, n=1):
    return (" " * n).join(t)


def make_wrap(fig):
    """Greedy word-wrap to a measured width, in figure fractions."""
    r = fig.canvas.get_renderer()

    def width(t, fs, style="normal"):
        p = fig.text(0, 0, t, fontsize=fs, style=style)
        w = p.get_window_extent(renderer=r).width / fig.bbox.width
        p.remove()
        return w

    def wrap(s, fs, max_w, style="normal"):
        out, line = [], ""
        for word in s.split():
            trial = word if not line else line + " " + word
            if line and width(trial, fs, style) > max_w:
                out.append(line)
                line = word
            else:
                line = trial
        if line:
            out.append(line)
        return "\n".join(out)

    return wrap


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    wrap = make_wrap(fig)
    y = 0.962

    # ---------------- header ----------------------------------------------
    ax.text(L, y, "How the sequence was interrogated", fontsize=17, color=INK,
            va="top", fontweight="bold")
    y -= lh(17) + 0.012
    ax.text(L, y,
            "Eight analyses on a ladder of abstraction. Each asks a different question, and "
            "each could have returned a positive result independently of the others.",
            fontsize=9.6, color=INK_2, va="top")
    y -= lh(9.6) + 0.016

    # column headers
    X_LEVEL = L
    X_TOOL, X_DESC, X_HIT = L + 0.150, L + 0.330, R - 0.215
    ax.text(X_TOOL, y, track("METHOD"), fontsize=7.0, color=INK_MUTED,
            va="top", fontweight="bold")
    ax.text(X_DESC, y, track("HOW IT WORKS"), fontsize=7.0, color=INK_MUTED,
            va="top", fontweight="bold")
    ax.text(X_HIT, y, track("IF POSITIVE, IT MEANT"), fontsize=7.0,
            color=INK_MUTED, va="top", fontweight="bold")
    y -= lh(7.0) + 0.008
    ax.plot([L, R], [y, y], color=INK_2, lw=1.0)
    y -= 0.016

    # ---------------- levels ----------------------------------------------
    ROW_GAP, LEVEL_GAP = 0.016, 0.026
    for li, (level, colour, tools) in enumerate(LEVELS):
        band_top = y
        ax.text(X_LEVEL, y, level, fontsize=8.2, color=colour, va="top",
                fontweight="bold", linespacing=1.38)
        for tool, how, hit in tools:
            # both prose columns wrap to their own measured width, so neither can
            # run into the next column when a column boundary moves
            w_how = wrap(how, 7.4, (X_HIT - X_DESC) - 0.022)
            w_hit = wrap(hit, 7.4, (R - X_HIT), style="italic")
            w_tool = wrap(tool, 8.6, (X_DESC - X_TOOL) - 0.018)
            ax.text(X_TOOL, y, w_tool, fontsize=8.6, color=INK, va="top",
                    fontweight="bold", linespacing=1.35)
            ax.text(X_DESC, y, w_how, fontsize=7.4, color=INK_2, va="top",
                    linespacing=1.42)
            ax.text(X_HIT, y, w_hit, fontsize=7.4, color=colour, va="top",
                    linespacing=1.42, style="italic")
            nlines = 1 + max(w_how.count("\n"), w_hit.count("\n"),
                             w_tool.count("\n"))
            y -= lh(7.4, 1.42) * (nlines - 1) + lh(7.4) + ROW_GAP
        # one accent rule spanning the level's rows, in its own gutter
        ax.add_patch(Rectangle((X_TOOL - 0.016, y + ROW_GAP - 0.002), 0.0035,
                               band_top - y - ROW_GAP + 0.002,
                               fc=colour, ec="none", zorder=3))
        y -= LEVEL_GAP - ROW_GAP
        if li < len(LEVELS) - 1:
            ax.plot([L, R], [y + LEVEL_GAP * 0.42, y + LEVEL_GAP * 0.42],
                    color=GRID, lw=0.8)

    # ---------------- the synthesis ---------------------------------------
    y -= 0.004
    pad = 0.018
    body = ("Every method above detects hazard, or origin, by resemblance to something already "
            "catalogued. They differ in what counts as resemblance — exact k-mers, scored "
            "substitutions, profile models, learned embeddings, predicted folds — and that is "
            "why running all eight is informative rather than redundant: a sequence can evade "
            "the strictest definition while tripping a looser one. A uniform negative across "
            "all eight levels is therefore a strong statement about the sequence's relationship "
            "to the known world, and no statement at all about what it does.")
    # wrap by measurement
    r = fig.canvas.get_renderer()

    def width(t, fs):
        p = fig.text(0, 0, t, fontsize=fs)
        w = p.get_window_extent(renderer=r).width / fig.bbox.width
        p.remove()
        return w

    maxw = (R - L) - 2 * pad
    lines, cur = [], ""
    for word in body.split():
        trial = word if not cur else cur + " " + word
        if cur and width(trial, 8.4) > maxw:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    lines.append(cur)
    wrapped = "\n".join(lines)
    panel_h = pad + lh(7.2) + 0.009 + lh(8.4, 1.5) * (len(lines) - 1) + lh(8.4) + pad
    ax.add_patch(FancyBboxPatch((L, y - panel_h), R - L, panel_h,
                                boxstyle="round,pad=0,rounding_size=0.008",
                                fc="#f3f4f6", ec="none", zorder=1))
    yy = y - pad
    ax.text(L + pad, yy, track("WHY EIGHT METHODS AND NOT ONE"), fontsize=7.2,
            color=INK_2, va="top", fontweight="bold", zorder=4)
    yy -= lh(7.2) + 0.009
    ax.text(L + pad, yy, wrapped, fontsize=8.4, color=INK, va="top",
            linespacing=1.5, zorder=4)
    y -= panel_h

    if y < 0.004:
        print(f"WARNING: overruns page bottom (cursor {y:+.4f})")
    else:
        print(f"bottom clear: cursor {y:+.4f}")

    fig.savefig(OUT, dpi=240, facecolor=SURFACE)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
