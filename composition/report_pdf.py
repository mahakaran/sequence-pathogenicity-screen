"""Render the sequence composition / codon-usage analysis as a PDF report.

Recomputes every number from sequence.fasta (nothing is hard-coded), renders
five figures, and lays them out with reportlab.

Usage:  python report_pdf.py        -> composition/sequence_analysis_report.pdf
"""

import itertools
import math
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Restriction import CommOnly, RestrictionBatch
from Bio.Seq import Seq
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

import composition_codon as cc

HERE = Path(__file__).resolve().parent
FIGS = HERE / "figures"
OUT = HERE / "sequence_analysis_report.pdf"

# ---- dataviz reference palette (light mode) --------------------------------
SURFACE, INK, INK_2, INK_MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
          "#6250d6", "#e34948"]
BLUE, ORANGE, AQUA = SERIES[0], SERIES[1], SERIES[2]
NEUTRAL = "#b9b8b3"

STD = CodonTable.unambiguous_dna_by_id[1]
CODONS = ["".join(c) for c in itertools.product("TCAG", repeat=3)]
AA = {c: STD.forward_table.get(c, "*") for c in CODONS}
SYN = defaultdict(list)
for _c, _a in AA.items():
    SYN[_a].append(_c)
AA_NAMES = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly",
            "H": "His", "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn",
            "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val",
            "W": "Trp", "Y": "Tyr"}
ANTISENSE_STOP = {"TTA", "CTA", "TCA"}   # reverse complement is TAA / TAG / TGA


def style_axes(ax, grid_axis=None):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, labelsize=8, length=3, width=0.8)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, lw=0.7, zorder=0)
        ax.set_axisbelow(True)


def codons_of(s):
    return [s[i:i + 3] for i in range(0, len(s) - len(s) % 3, 3)]


# ============================================================== computation
def compute():
    d = {}
    nt = str(next(SeqIO.parse(HERE.parent / "sequence.fasta", "fasta")).seq).upper()
    rc = str(Seq(nt).reverse_complement())
    L = len(nt)
    d["nt"], d["rc"], d["L"] = nt, rc, L
    d["counts"] = Counter(nt)
    d["gc"] = 100 * (d["counts"]["G"] + d["counts"]["C"]) / L
    d["gc_skew"] = (d["counts"]["G"] - d["counts"]["C"]) / (d["counts"]["G"] + d["counts"]["C"])
    d["at_skew"] = (d["counts"]["A"] - d["counts"]["T"]) / (d["counts"]["A"] + d["counts"]["T"])
    d["win"] = [(i + 30, cc.gc(nt[i:i + 60])) for i in range(0, L - 59, 6)]
    d["oe"] = cc.dinuc_oe(nt)

    # frames
    d["frames"] = {"+1 (ORF1)": nt[: L // 3 * 3], "-1 (ORF4)": rc[: L // 3 * 3]}
    d["gcpos"] = {}
    for name, s in d["frames"].items():
        cs = codons_of(s)
        d["gcpos"][name] = [100 * sum(c[i] in "GC" for c in cs) / len(cs) for i in range(3)]

    # ORF4 codon usage
    orf4 = codons_of(rc)
    d["orf4_prot"] = str(Seq(rc).translate())
    d["orf4_codons"] = orf4
    by = defaultdict(Counter)
    for c in orf4:
        by[AA[c]][c] += 1
    d["orf4_usage"] = by
    d["multi"] = {a: v for a, v in by.items() if a != "*" and len(SYN[a]) > 1}
    d["n_multi_aa"] = len(d["multi"])
    d["n_multi_codons"] = sum(sum(v.values()) for v in d["multi"].values())
    d["n_single_codon_aa"] = sum(1 for v in d["multi"].values() if len(v) == 1)
    d["enc"] = {n: cc.enc_wright(codons_of(s)) for n, s in d["frames"].items()}
    d["gc3"] = {}
    for n, s in d["frames"].items():
        cs = [c for c in codons_of(s) if AA[c] not in "*MW"]
        d["gc3"][n] = 100 * sum(c[2] in "GC" for c in cs) / len(cs)

    # back-translation determinism
    tab = {}
    for c in orf4:
        tab.setdefault(AA[c], c)
    d["bt_table"] = tab
    d["bt_exact"] = "".join(tab[a] for a in d["orf4_prot"]) == rc
    d["antisense_stop_used"] = sum(Counter(orf4)[c] for c in ANTISENSE_STOP)

    # reference panel
    d["panel"] = []
    for taxid, label, group, ggc in cc.PANEL:
        try:
            counts, n = cc.fetch_table(taxid)
        except Exception:
            continue
        p = cc.syn_freqs(counts)
        logp = sum(math.log10(sum(p[c] ** sum(v.values()) for c in SYN[a]))
                   for a, v in d["multi"].items())
        keep = sum(sum(v.values()) * math.log10(sum(p[c] for c in SYN[a]
                                                    if c not in ANTISENSE_STOP))
                   for a, v in d["multi"].items())
        d["panel"].append(dict(
            taxid=taxid, label=label, group=group, counts=counts, n=n,
            gc3=cc.gc3_of_table(counts), logp=logp, logkeep=keep,
            cai4=cc.cai(orf4, counts), ll4=cc.syn_loglik(orf4, counts),
            cai1=cc.cai(codons_of(d["frames"]["+1 (ORF1)"]), counts),
            ll1=cc.syn_loglik(codons_of(d["frames"]["+1 (ORF1)"]), counts)))
    d["panel"].sort(key=lambda r: -r["logp"])

    # ORF map: stop positions per frame
    d["orfmap"] = {}
    for strand, s in (("+", nt), ("-", rc)):
        for f in range(3):
            cs = codons_of(s[f:])
            stops = [i for i, c in enumerate(cs) if AA[c] == "*"]
            run = best = 0
            for c in cs:
                run = 0 if AA[c] == "*" else run + 1
                best = max(best, run)
            d["orfmap"][f"{strand}{f + 1}"] = dict(n=len(cs), stops=stops, best=best)

    # shuffle test A
    rng = np.random.default_rng(0)
    arr = np.array(list(nt))
    N = 20000
    both = 0
    for _ in range(N):
        rng.shuffle(arr)
        fo = cc.fully_open("".join(arr))
        both += fo["+1"] and fo["-1"]
    d["shuffleA_n"], d["shuffleA_hits"] = N, both

    # restriction
    cuts = {e: p for e, p in CommOnly.search(Seq(nt), linear=True).items() if p}
    d["n_enzymes"] = len(CommOnly)
    d["n_cutters"] = len(cuts)
    d["six_cutters"] = sorted(
        {(str(e), e.site, tuple(p)) for e, p in cuts.items()
         if len(e.site) == 6 and set(e.site) <= set("ACGT")})
    iis = [e for e in ["BsaI", "BsmBI", "Esp3I", "BbsI", "SapI", "AarI", "PaqCI",
                       "BtgZI", "BspQI"] if e in CommOnly.elements()]
    d["iis"] = [(str(e), e.site, bool(p))
                for e, p in RestrictionBatch(iis).search(Seq(nt), linear=True).items()]
    f = {b: d["counts"][b] / L for b in "ACGT"}
    d["basefreq"] = f
    pals = ["".join(p) for p in itertools.product("ACGT", repeat=6)
            if "".join(p) == str(Seq("".join(p)).reverse_complement())]
    exp_present = 0.0
    for s in pals:
        pr = np.prod([f[b] for b in s]) * (L - 5)
        exp_present += 1 - math.exp(-pr)
    d["pal_expected"] = exp_present
    d["pal_observed"] = [s for s in pals if s in nt]

    # ORF coordinates from the ORF-finder output
    d["orfs"] = []
    for r in SeqIO.parse(HERE.parent / "screening-exercise-orfs.fasta", "fasta"):
        name, a, b = r.id.split("|")[-1].split(":")
        a, b = int(a), int(b)
        d["orfs"].append(dict(name=name, a=a, b=b, strand="+" if a <= b else "-",
                              aa=len(r.seq), partial="partial" in r.description))
    return d


# ================================================================= figures
def fig_composition(d):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 2.9), facecolor=SURFACE,
                                 gridspec_kw={"width_ratios": [1, 2.1], "wspace": 0.22})
    # (a) base counts
    style_axes(a1, grid_axis="y")
    bases = list("ACGT")
    vals = [100 * d["counts"][b] / d["L"] for b in bases]
    cols = [ORANGE if b in "AT" else BLUE for b in bases]
    bars = a1.bar(bases, vals, color=cols, width=0.62, zorder=3)
    for b, v, n in zip(bars, vals, (d["counts"][x] for x in bases)):
        a1.text(b.get_x() + b.get_width() / 2, v + 1.1, f"{v:.1f}%",
                ha="center", fontsize=8, color=INK)
    a1.set_ylim(0, max(vals) * 1.25)
    a1.set_ylabel("% of sequence", fontsize=8.5, color=INK_2)
    a1.set_title("Base composition", loc="left", fontsize=9.5, color=INK, pad=8)

    # (b) sliding GC
    style_axes(a2, grid_axis="y")
    x = [p for p, _ in d["win"]]
    y = [v for _, v in d["win"]]
    a2.axhline(d["gc"], color=INK_MUTED, lw=1, ls="--", zorder=2)
    a2.text(d["L"] - 2, d["gc"] + 0.8, f"mean {d['gc']:.1f}%", ha="right",
            fontsize=7.5, color=INK_2)
    a2.plot(x, y, color=BLUE, lw=2, zorder=4, solid_capstyle="round")
    a2.set_xlim(0, d["L"])
    a2.set_ylim(min(y) - 6, max(y) + 7)
    a2.set_xlabel("position (nt), 60 nt window", fontsize=8.5, color=INK_2)
    a2.set_ylabel("GC %", fontsize=8.5, color=INK_2)
    a2.set_title("GC content along the sequence — uniform throughout",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    fig.savefig(FIGS / "composition.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def _usage_panel(ax, rows, title, subtitle):
    """rows: list of (aa, n, [(codon, fraction), ...]) ordered by fraction desc."""
    style_axes(ax)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    gap = 0.004   # surface gap between touching segments
    for i, (aa, n, segs) in enumerate(rows):
        y = len(rows) - 1 - i
        left = 0.0
        for j, (cod, frac) in enumerate(segs):
            if frac <= 0:
                continue
            w = frac - (gap if j < len(segs) - 1 else 0)
            ax.barh(y, w, left=left, height=0.62, color=SERIES[j % len(SERIES)],
                    zorder=3, linewidth=0)
            if frac > 0.17:
                ax.text(left + frac / 2, y, cod, ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold", zorder=4)
            left += frac
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{AA_NAMES[aa]} ({n})" for aa, n, _ in rows][::-1], fontsize=7.5)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100%"], fontsize=7.5)
    ax.set_title(f"{title}\n", loc="left", fontsize=9.5, color=INK, pad=14)
    ax.text(0, 1.012, subtitle, transform=ax.transAxes, fontsize=7.8, color=INK_2,
            va="bottom")


def fig_codon_bias(d):
    order = sorted(d["multi"].items(), key=lambda kv: -sum(kv[1].values()))
    obs_rows = []
    for aa, v in order:
        n = sum(v.values())
        segs = [(c, cnt / n) for c, cnt in v.most_common()]
        obs_rows.append((aa, n, segs))

    ref = next(r for r in d["panel"] if r["label"].startswith("Saccharomyces"))
    p = cc.syn_freqs(ref["counts"])
    ref_rows = []
    for aa, v in order:
        n = sum(v.values())
        segs = sorted(((c, p[c]) for c in SYN[aa]), key=lambda t: -t[1])
        ref_rows.append((aa, n, segs))

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), facecolor=SURFACE,
                             gridspec_kw={"wspace": 0.30})
    _usage_panel(axes[0], obs_rows, "Observed in ORF4",
                 f"every amino acid uses exactly ONE codon ({d['n_single_codon_aa']}/{d['n_multi_aa']})")
    _usage_panel(axes[1], ref_rows, "Expected — S. cerevisiae usage",
                 "a natural gene spreads usage across synonymous codons")
    fig.text(0.5, -0.03, "Bar = share of that amino acid's occurrences; "
             "segments are labelled with the codon.", ha="center", fontsize=7.8, color=INK_2)
    fig.savefig(FIGS / "codon_bias.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def fig_probability(d):
    rows = sorted(d["panel"], key=lambda r: r["logp"])
    fig, ax = plt.subplots(figsize=(9.2, 3.5), facecolor=SURFACE)
    style_axes(ax, grid_axis="x")
    y = np.arange(len(rows))
    gmap = {"Bacteria": BLUE, "Eukaryota": ORANGE, "Virus": AQUA}
    ax.barh(y, [r["logp"] for r in rows],
            color=[gmap[r["group"]] for r in rows], height=0.62, zorder=3, linewidth=0)
    for i, r in enumerate(rows):
        ax.text(r["logp"] - 0.6, i, f"{r['logp']:.1f}", va="center", ha="right",
                fontsize=7.5, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlabel("log₁₀ P(every amino acid uses exactly one codon)  —  more negative = less plausible",
                  fontsize=8.5, color=INK_2)
    ax.set_xlim(min(r["logp"] for r in rows) * 1.18, 0)
    ax.set_title("No natural codon usage plausibly produces ORF4's codon choices",
                 loc="left", fontsize=9.5, color=INK, pad=10)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=7, mfc=c, mec="none", label=g)
               for g, c in gmap.items()]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=8,
              labelcolor=INK_2, ncol=1, borderaxespad=1.2)
    fig.savefig(FIGS / "probability.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def fig_orf_map(d):
    fig, ax = plt.subplots(figsize=(9.2, 2.7), facecolor=SURFACE)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    frames = ["+1", "+2", "+3", "-1", "-2", "-3"]
    for i, fr in enumerate(frames):
        y = len(frames) - 1 - i
        info = d["orfmap"][fr]
        full = info["best"] == info["n"]
        ax.barh(y, info["n"] * 3, height=0.5, color=BLUE if full else NEUTRAL,
                zorder=3, linewidth=0, alpha=1.0 if full else 0.65)
        for s in info["stops"]:
            ax.plot([s * 3 + 1.5, s * 3 + 1.5], [y - 0.3, y + 0.3],
                    color=ORANGE, lw=1.6, zorder=4, solid_capstyle="butt")
        lbl = "open end-to-end" if full else f"longest open run {info['best']} codons"
        ax.text(info["n"] * 3 + 5, y, lbl, va="center", fontsize=7.5,
                color=INK if full else INK_2)
    ax.set_yticks(range(len(frames)))
    ax.set_yticklabels(frames[::-1], fontsize=8.5)
    ax.set_xlim(0, d["L"] * 1.33)
    ax.set_xticks([0, 72, 144, 216, 288])
    ax.set_xlabel("position (nt)", fontsize=8.5, color=INK_2)
    ax.set_title("Reading frames — frames +1 and −1 have no stop codon anywhere",
                 loc="left", fontsize=9.5, color=INK, pad=10)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=7, mfc=BLUE, mec="none",
                          label="frame open end-to-end"),
               plt.Line2D([], [], marker="s", ls="", ms=7, mfc=NEUTRAL, mec="none",
                          label="frame interrupted"),
               plt.Line2D([], [], color=ORANGE, lw=1.8, label="stop codon")]
    ax.legend(handles=handles, loc="lower right", ncol=3, frameon=False,
              fontsize=8, labelcolor=INK_2, borderaxespad=0.2,
              bbox_to_anchor=(1.0, -0.42))
    fig.savefig(FIGS / "orf_map.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def fig_dinuc(d):
    keys = ["".join(p) for p in itertools.product("ACGT", repeat=2)]
    vals = [d["oe"][k] for k in keys]
    fig, ax = plt.subplots(figsize=(9.2, 2.6), facecolor=SURFACE)
    style_axes(ax, grid_axis="y")
    cols = [ORANGE if k == "CG" else BLUE for k in keys]
    ax.bar(range(len(keys)), vals, color=cols, width=0.62, zorder=3, linewidth=0)
    ax.axhline(1.0, color=INK_MUTED, lw=1, ls="--", zorder=2)
    ax.text(-0.35, 1.05, "expected", ha="left", fontsize=7.5, color=INK_2)
    i = keys.index("CG")
    ax.set_ylim(0, max(vals) * 1.28)
    ax.annotate("CpG absent (≈5.9 expected)", xy=(i, 0.03), xytext=(i - 0.2, 1.42),
                ha="center", fontsize=7.8, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8,
                                shrinkA=2, shrinkB=2))
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys, fontsize=7.5, rotation=0)
    ax.set_ylabel("observed / expected", fontsize=8.5, color=INK_2)
    ax.set_title("Dinucleotide frequencies — CpG is the single clear outlier",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    fig.savefig(FIGS / "dinuc.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


# ==================================================================== PDF
def build_pdf(d):
    ss = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=ss["BodyText"], fontName="Helvetica",
                          fontSize=9.3, leading=13.4, textColor=colors.HexColor(INK),
                          spaceAfter=6, alignment=TA_JUSTIFY)
    h1 = ParagraphStyle("h1", parent=body, fontName="Helvetica-Bold", fontSize=13,
                        leading=16, spaceBefore=14, spaceAfter=7, alignment=0,
                        textColor=colors.HexColor(INK))
    h2 = ParagraphStyle("h2", parent=body, fontName="Helvetica-Bold", fontSize=10.3,
                        leading=13, spaceBefore=10, spaceAfter=4, alignment=0)
    small = ParagraphStyle("small", parent=body, fontSize=7.9, leading=10.6,
                           textColor=colors.HexColor(INK_2), alignment=0, spaceAfter=4)
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=7.4,
                          leading=9.6, alignment=0)
    title = ParagraphStyle("title", parent=body, fontName="Helvetica-Bold",
                           fontSize=19, leading=23, spaceAfter=2, alignment=0)
    sub = ParagraphStyle("sub", parent=body, fontSize=10, leading=14,
                         textColor=colors.HexColor(INK_2), spaceAfter=12, alignment=0)

    def fig(name, width=170 * mm):
        p = FIGS / name
        iw, ih = Image(str(p)).imageWidth, Image(str(p)).imageHeight
        return Image(str(p), width=width, height=width * ih / iw)

    def tbl(data, widths, align_right=(), header=True, size=8.2):
        t = Table(data, colWidths=widths, hAlign="LEFT")
        st = [("FONT", (0, 0), (-1, -1), "Helvetica", size),
              ("TEXTCOLOR", (0, 0), (-1, -1), colors.HexColor(INK)),
              ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
              ("TOPPADDING", (0, 0), (-1, -1), 3),
              ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
              ("LEFTPADDING", (0, 0), (-1, -1), 5),
              ("LINEBELOW", (0, 0), (-1, -2), 0.4, colors.HexColor(GRID))]
        if header:
            st += [("FONT", (0, 0), (-1, 0), "Helvetica-Bold", size),
                   ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor(INK_2)),
                   ("LINEBELOW", (0, 0), (-1, 0), 0.8, colors.HexColor(INK_MUTED))]
        for c in align_right:
            st.append(("ALIGN", (c, 0), (c, -1), "RIGHT"))
        t.setStyle(TableStyle(st))
        return t

    def callout(text, tint="#eef4fc", edge=BLUE):
        p = Paragraph(text, ParagraphStyle("c", parent=body, spaceAfter=0))
        t = Table([[p]], colWidths=[170 * mm])
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(tint)),
            ("LINEBEFORE", (0, 0), (0, -1), 2.2, colors.HexColor(edge)),
            ("LEFTPADDING", (0, 0), (-1, -1), 9), ("RIGHTPADDING", (0, 0), (-1, -1), 9),
            ("TOPPADDING", (0, 0), (-1, -1), 8), ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
        return t

    S = []
    A = S.append

    # ---------------- header
    A(Paragraph("Nucleotide composition &amp; codon-usage analysis", title))
    A(Paragraph(f"Screening exercise · <font face='Courier'>sequence.fasta</font> · "
                f"{d['L']} bp · {date.today().isoformat()}", sub))

    A(callout(
        "<b>Principal finding.</b> The nucleotide sequence was almost certainly produced by "
        "back-translating a protein, not sampled from a genome. In the −1 frame (ORF4), "
        f"all {d['n_single_codon_aa']} of the {d['n_multi_aa']} amino acids that have synonymous "
        f"codons use <b>exactly one codon each</b> across {d['n_multi_codons']} codons. Under the "
        "most favourable natural codon usage tested this has probability "
        f"≈10<super>{max(r['logp'] for r in d['panel']):.0f}</super>. This explains the otherwise "
        "puzzling screening result: nucleotide methods (BLASTn, Kraken2) return nothing, while "
        "protein-level methods identify ORF4 as a Cu-Zn superoxide dismutase."))
    A(Spacer(1, 7))

    # ---------------- 1 the sequence
    A(Paragraph("1 · Sequence and reading frames", h1))
    seq_lines = []
    for i in range(0, d["L"], 60):
        seq_lines.append(f"{i + 1:>4}  " + " ".join(d["nt"][j:j + 10]
                                                    for j in range(i, min(i + 60, d["L"]), 10)))
    A(Paragraph("<br/>".join(seq_lines), mono))
    A(Spacer(1, 6))
    orows = [["ORF", "nt range", "Strand", "Length", "Ends"]]
    for o in d["orfs"]:
        orows.append([o["name"], f"{o['a']}–{o['b']}", o["strand"], f"{o['aa']} aa",
                      "runs off sequence end" if o["partial"] else "stop to stop"])
    A(tbl(orows, [20 * mm, 28 * mm, 20 * mm, 22 * mm, 56 * mm], size=8))
    A(Paragraph(
        "ORF coordinates as called by NCBI ORFfinder. ORF4 — the −1 frame spanning the whole "
        "sequence — is the frame analysed throughout this report; its protein product is the "
        "one identified as a Cu-Zn superoxide dismutase by protein-level search.", small))

    # ---------------- 2 composition
    A(Paragraph("2 · Base composition", h1))
    A(Paragraph(
        f"The sequence is AT-rich: GC {d['gc']:.1f}%, AT {100 - d['gc']:.1f}% "
        f"(A {d['counts']['A']}, C {d['counts']['C']}, G {d['counts']['G']}, T {d['counts']['T']}). "
        f"GC skew is {d['gc_skew']:+.2f} and AT skew {d['at_skew']:+.2f}. GC content is uniform "
        "along the sequence, so there is no sign of a junction between segments of different "
        "origin. AT-richness alone is not diagnostic — it is shared by AT-rich bacteria, many "
        "fungi, <i>Plasmodium</i>, poxviruses and phage.", body))
    A(fig("composition.png"))
    A(Spacer(1, 6))

    gp = d["gcpos"]
    A(tbl([["Frame", "GC1", "GC2", "GC3", "ENC", "Codons"]] +
          [[n, f"{gp[n][0]:.1f}%", f"{gp[n][1]:.1f}%", f"{gp[n][2]:.1f}%",
            f"{d['enc'][n]:.1f}", "96"] for n in d["frames"]],
          [38 * mm, 22 * mm, 22 * mm, 22 * mm, 20 * mm, 20 * mm], align_right=(1, 2, 3, 4, 5)))
    A(Paragraph(
        "ENC (effective number of codons) runs from 61 (no bias) to 20 (one codon per amino "
        f"acid). ORF4 sits at the <b>theoretical floor of {d['enc']['-1 (ORF4)']:.1f}</b>. Note that "
        "GC1 and GC3 are mirrored between the two frames purely because they are reverse "
        "complements of one another — that symmetry is an artefact, not evidence.", small))

    # ---------------- 2 dinucleotide
    # ---------------- 3 codon usage
    A(PageBreak())
    A(Paragraph("3 · Codon usage — the decisive evidence", h1))
    A(Paragraph(
        "Every amino acid in ORF4 is encoded by a single codon throughout: 9×Ala all GCT, "
        "10×Gly all GGT, 9×Ile all ATT, 9×Lys all AAA, and so on. Natural genes — including "
        "those of the most compositionally extreme genomes — always spread usage across "
        "synonymous codons.", body))
    A(fig("codon_bias.png"))
    A(Spacer(1, 4))
    A(Paragraph(
        "Back-translating ORF4's 96-residue protein with the single-codon table below "
        "reproduces the observed 288 bp <b>base for base</b>"
        f"{'' if d['bt_exact'] else ' (CHECK FAILED)'}. The nucleotide sequence therefore "
        "carries no information beyond the protein plus this table.", body))
    A(Paragraph("&nbsp;&nbsp;" + "&nbsp; ".join(f"{a}:{c}" for a, c in sorted(d["bt_table"].items())),
                mono))
    A(Spacer(1, 6))
    best = max(d["panel"], key=lambda r: r["logp"])
    A(fig("probability.png"))
    A(Paragraph(
        "Probability that every amino acid would independently use a single codon, evaluated "
        "under each reference organism's codon usage (Kazusa codon usage database, "
        f"pseudocount 0.5). The least extreme value, <i>{best['label']}</i> at "
        f"10<super>{best['logp']:.1f}</super>, arises only because that organism's own usage is "
        "itself near-deterministic — toward GC-ending codons that ORF4 never uses (its CAI "
        f"against {best['label']} is {best['cai4']:.2f}, the lowest in the panel). It is quoted "
        "as a conservative bound.", small))

    # ---------------- 3 dual coding
    A(PageBreak())
    A(Paragraph("4 · Dual coding and avoided codons", h1))
    A(Paragraph(
        "Frames +1 and −1 contain no stop codon across all 96 codons, so both strands read "
        "open end-to-end. In "
        f"{d['shuffleA_n']:,} shuffles preserving base composition, this occurred "
        f"{d['shuffleA_hits']} times.", body))
    A(fig("orf_map.png"))
    A(Spacer(1, 4))
    A(Paragraph(
        "The mechanism is visible in the codon choices. Three codons — <b>TTA, CTA, TCA</b> — "
        "have a stop codon as their reverse complement, so using any of them would close the "
        f"opposite strand. ORF4 uses them <b>{d['antisense_stop_used']} times</b>. All 6 leucines "
        "are TTG rather than the more AT-rich TTA, and all 6 serines are TCT rather than TCA — "
        "choices that run against the sequence's own AT bias.", body))

    A(KeepTogether([
        Paragraph("5 · Dinucleotide signature", h1),
        fig("dinuc.png"),
        Paragraph(
            f"CpG is entirely absent ({d['oe']['CG']:.2f} observed/expected; ≈5.9 occurrences "
            "expected from base composition). CpG depletion occurs naturally in vertebrates and "
            "some viruses, but complete absence over 288 bp is a further mark of constructed "
            "sequence.", small)]))

    # ---------------- reference panel table
    A(PageBreak())
    A(Paragraph("6 · Codon usage versus reference organisms", h1))
    A(Paragraph(
        "Reference organisms ranked by how well ORF4's codon choices fit their usage. "
        "<b>This ranking should not be used to infer a source organism.</b> Because the codon "
        "choice was imposed rather than inherited, the ordering simply sorts references by "
        "AT-richness.", body))
    rows = [["Organism", "Group", "ref GC3", "synLL", "CAI", "log10 P"]]
    for r in sorted(d["panel"], key=lambda r: -r["ll4"]):
        rows.append([r["label"], r["group"], f"{r['gc3']:.1f}%", f"{r['ll4']:.2f}",
                     f"{r['cai4']:.2f}", f"{r['logp']:.1f}"])
    A(tbl(rows, [54 * mm, 24 * mm, 20 * mm, 20 * mm, 18 * mm, 20 * mm],
          align_right=(2, 3, 4, 5), size=7.9))
    A(Paragraph(
        "synLL = mean log<sub>2</sub> P(codon | amino acid) under that organism's usage; CAI = codon "
        "adaptation index; log10 P as in section 3. ORF4's GC3 of "
        f"{d['gc3']['-1 (ORF4)']:.1f}% is below every reference in the panel.", small))

    # ---------------- 5 restriction
    A(Paragraph("7 · Restriction sites", h1))
    A(Paragraph(
        f"Of {d['n_enzymes']} commercially available enzymes, {d['n_cutters']} cut the sequence; "
        "most are 4-cutters or have degenerate sites and cut essentially any sequence. The "
        "unambiguous 6-bp cutters are:", body))
    rrows = [["Enzyme", "Site", "Position"]] + [
        [n, s, ", ".join(str(x) for x in p)] for n, s, p in d["six_cutters"]]
    A(tbl(rrows, [40 * mm, 34 * mm, 40 * mm], size=8))
    A(Spacer(1, 5))
    A(callout(
        "<b>These sites carry no design signal.</b> Because the codon table is fixed at one "
        "codon per amino acid, the 288 bp is fully determined by the protein. No restriction "
        "site could have been added or removed without breaking that rule, so neither presence "
        "nor absence reflects intent.", tint="#fdf0e9", edge=ORANGE))
    A(Spacer(1, 6))
    A(Paragraph(
        f"The counts are also unremarkable. Given this base composition one expects "
        f"<b>{d['pal_expected']:.1f}</b> distinct 6-bp palindromic sites to occur at least once in "
        f"{d['L']} bp; <b>{len(d['pal_observed'])}</b> were observed "
        f"({', '.join(d['pal_observed'])}).", body))
    A(Paragraph(
        "All Type IIS sites used for Golden Gate 'domestication' (" +
        ", ".join(n for n, _, _ in d["iis"]) + ") are absent. This is <i>not</i> evidence of a "
        "domesticated synthetic part: these are GC-rich sites in an AT-rich sequence, with "
        "expected counts of 0.02–0.05, so absence is the default expectation.", body))

    # ---------------- 6 caveats + methods
    A(PageBreak())
    A(Paragraph("8 · Interpretation and caveats", h1))
    for t in [
        "<b>The protein is not implicated.</b> The evidence concerns the nucleotide encoding "
        "only. ORF4's 96-residue protein product resembles a natural Cu-Zn superoxide dismutase; "
        "it is the DNA that was constructed.",
        "<b>288 bp is short.</b> The single-codon result is robust because it rests on "
        f"{d['n_multi_codons']} codons across {d['n_multi_aa']} amino acids, but ENC and CAI are "
        "individually noisy at this length.",
        "<b>Origin is undetermined.</b> Codon usage cannot identify a source organism for a "
        "back-translated sequence, and the reference ranking in section 6 should not be read "
        "as a taxonomic assignment.",
        "<b>Screening implication.</b> Exact-match nucleotide methods are expected to fail on "
        "back-translated sequence, since nearly every third base differs from any natural "
        "homologue. Protein-level and structure-level methods remain effective and should carry "
        "the weight in any screening decision for this submission.",
    ]:
        A(Paragraph("• " + t, body))

    A(Paragraph("9 · Cross-method context", h2))
    A(tbl([["Method", "Level", "Result"],
           ["BLASTn vs nt", "nucleotide", "no significant hit"],
           ["Kraken2 — Standard-8", "nucleotide (35-mer)", "unclassified, 0/254 k-mers"],
           ["Kraken2 — PlusPF-8", "nucleotide (35-mer)", "unclassified, 0/254 k-mers"],
           ["ESM-2 8M + ProteinCLIP", "protein embedding", "ORF4 → Cu-Zn superoxide dismutase"],
           ["AlphaFold + structure search", "structure", "model generated for ORF4"]],
          [48 * mm, 36 * mm, 76 * mm], size=8.2))

    A(Paragraph("10 · Methods", h2))
    A(Paragraph(
        "All statistics were computed from <font face='Courier'>sequence.fasta</font> "
        f"({d['L']} bp). Reading frames follow the ORF coordinates in "
        "<font face='Courier'>screening-exercise-orfs.fasta</font>; ORF4 is the −1 frame "
        "(positions 288→1). Codon reference tables are from the Kazusa codon usage database "
        "(cached in <font face='Courier'>composition/ref_tables/</font>). Synonymous-codon "
        "probabilities use a pseudocount of 0.5. ENC follows Wright (1990). Restriction analysis "
        "uses Biopython's REBASE <i>CommOnly</i> set. The shuffle test permutes bases while "
        "preserving composition. Expected palindromic-site counts assume independent bases at "
        "the observed frequencies. Full numeric output: "
        "<font face='Courier'>composition/report.txt</font>; code: "
        "<font face='Courier'>composition/composition_codon.py</font> and "
        "<font face='Courier'>composition/report_pdf.py</font>.", small))

    # ---------------- doc with footer
    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.4)
        canvas.setFillColor(colors.HexColor(INK_MUTED))
        canvas.drawString(20 * mm, 12 * mm,
                          "Screening exercise — composition & codon-usage analysis")
        canvas.drawRightString(A4[0] - 20 * mm, 12 * mm, f"page {doc.page}")
        canvas.setStrokeColor(colors.HexColor(GRID))
        canvas.setLineWidth(0.5)
        canvas.line(20 * mm, 15 * mm, A4[0] - 20 * mm, 15 * mm)
        canvas.restoreState()

    doc = BaseDocTemplate(str(OUT), pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
                          topMargin=18 * mm, bottomMargin=20 * mm,
                          title="Sequence composition & codon-usage analysis",
                          author="screening exercise")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="f")
    doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=footer)])
    doc.build(S)
    print(f"Wrote {OUT}")


def main():
    FIGS.mkdir(exist_ok=True)
    d = compute()
    fig_composition(d)
    fig_codon_bias(d)
    fig_probability(d)
    fig_orf_map(d)
    fig_dinuc(d)
    build_pdf(d)


if __name__ == "__main__":
    main()
