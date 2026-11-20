"""Top-50 nearest neighbours per ORF, as table figures.

Ranked by cosine distance in the 1280-d ESM-2 650M embedding space (NOT by UMAP
distance, which distorts local geometry). Each row also carries the DIAMOND
e-value for that same pair, so embedding proximity and actual sequence homology
can be compared directly.

Pathogenicity is reported for the NEIGHBOUR protein, never for the ORF:
  VF      homology to a VFDB curated virulence factor (DIAMOND, e<=1e-5, id>=30%)
  Path    source organism is a recognised human/animal pathogen

Usage:  python neighbour_tables.py
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
RES = HERE / "results"
DIA = HERE / "diamond"
FIGDIR = RES / "neighbour_tables"

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
ORF_COLOR = "#eb6834"
HIT = "#e34948"        # pathogen / virulence flag
OK = "#1baf7a"

# Genera whose members are recognised human or animal pathogens. Hand-curated for the
# 54 genera that actually appear among the neighbours, cross-checked against the
# organism strings in Swiss-Prot (many carry explicit pathovar/serotype names).
PATHOGEN_GENERA = {
    "Acinetobacter", "Bacillus", "Bacteroides", "Campylobacter", "Chlamydia",
    "Citrobacter", "Clostridium", "Escherichia", "Klebsiella", "Mycoplasma",
    "Pseudomonas", "Rickettsia", "Shigella", "Streptococcus", "Yersinia",
}
# Venomous animals: their toxins are hazardous but they are not infectious agents,
# so they are flagged separately rather than as pathogens.
VENOMOUS_GENERA = {"Bitis", "Crotalus", "Micrurus", "Naja", "Tityus", "Lychas"}
# Viral organism names (not binomial) matched by substring.
PATHOGEN_VIRUS = [
    "variola", "monkeypox", "vaccinia", "varicella", "herpesvirus", "influenza",
    "rotavirus", "adenovirus", "papillomavirus", "coronavirus", "mimivirus",
]


def abbreviate(name, limit=46):
    """Shorten a protein name for the table without losing its distinguishing part."""
    if not isinstance(name, str):
        return "—"
    s = re.sub(r"\s*\((Fragment|Precursor)\)\s*", "", name).strip()
    for long, short in [("Putative uncharacterized protein", "Putative uncharacterized"),
                        ("Uncharacterized protein", "Uncharacterized"),
                        ("superfamily", "sf."), ("domain-containing", "domain-cont.")]:
        s = s.replace(long, short)
    return s if len(s) <= limit else s[:limit - 1] + "…"


def pathogen_status(organism, dom):
    if not isinstance(organism, str):
        return ""
    low = organism.lower()
    if dom == "Viruses":
        return "yes" if any(v in low for v in PATHOGEN_VIRUS) else ""
    genus = organism.split()[0]
    if genus in PATHOGEN_GENERA:
        return "yes"
    if genus in VENOMOUS_GENERA:
        return "venom"
    return ""


def build():
    d = pd.read_parquet(RES / "top50_neighbours.parquet")
    sys.path.insert(0, str(HERE))
    from annotate import domain, family

    meta = pq.read_table(PROJECT / "DS569k.parquet",
                         columns=["accession", "protein_name", "organism_name",
                                  "sequence_length", "ncbi_taxonomy_phylum",
                                  "ncbi_taxonomy_class"]).to_pandas().set_index("accession")
    coords = pd.read_parquet(RES / "umap650_coords.parquet").set_index("accession")

    j = meta.reindex(d.accession)
    d["protein_name"] = j.protein_name.values
    d["organism"] = j.organism_name.values
    d["length"] = j.sequence_length.values
    d["domain"] = [domain(p, c, o) for p, c, o in
                   zip(j.ncbi_taxonomy_phylum, j.ncbi_taxonomy_class, j.organism_name)]
    d["superfamily"] = coords.reindex(d.accession).superfamily.values
    d["family"] = [family(n) if isinstance(n, str) else None for n in d.protein_name]

    # DIAMOND e-value for the ORF <-> neighbour pair (empty file = no alignment at all)
    ev = {}
    p = DIA / "orf_vs_neighbours.tsv"
    if p.exists() and p.stat().st_size:
        t = pd.read_csv(p, sep="\t", header=None,
                        names=["q", "s", "pident", "length", "evalue", "bits"])
        for _, r in t.iterrows():
            key = (r.q.split("|")[-1].split(":")[0], r.s)
            ev[key] = min(ev.get(key, np.inf), r.evalue)
    d["diamond_e"] = [ev.get((o, a)) for o, a in zip(d.orf, d.accession)]

    # VFDB virulence-factor homology for the NEIGHBOUR protein
    vf = {}
    p = DIA / "neighbours_vs_vfdb.tsv"
    if p.exists() and p.stat().st_size:
        t = pd.read_csv(p, sep="\t", header=None,
                        names=["q", "s", "pident", "evalue", "bits", "stitle"])
        for _, r in t.iterrows():
            cat = re.search(r"-\s*([^()\[\]]+?)\s*\(VFC\d+\)", str(r.stitle))
            vf[r.q] = (cat.group(1).strip() if cat else "virulence factor", r.evalue)
    d["vf"] = [vf.get(a, ("", None))[0] for a in d.accession]
    d["vf_e"] = [vf.get(a, ("", None))[1] for a in d.accession]
    d["pathogen"] = [pathogen_status(o, dm) for o, dm in zip(d.organism, d.domain)]
    d.to_parquet(RES / "top50_annotated.parquet", index=False)
    return d


def render(d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGDIR.mkdir(parents=True, exist_ok=True)
    orfs = list(dict.fromkeys(d.orf))
    seqs = {}
    for line in (PROJECT / "screening-exercise-orfs.fasta").read_text().split(">"):
        if line.strip():
            h, *rest = line.strip().split("\n")
            seqs[h.split("|")[-1].split(":")[0]] = "".join(rest)

    cols = [("#", 0.030), ("Accession", 0.070), ("Cos.dist", 0.056), ("DIAMOND e", 0.068),
            ("Protein (abbreviated)", 0.270), ("Superfamily", 0.185),
            ("Organism", 0.165), ("Dom", 0.052), ("Path", 0.044), ("VF", 0.060)]

    for orf in orfs:
        sub = d[d.orf == orf].sort_values("rank")
        fig_h = 1.9 + 0.213 * len(sub)
        fig, ax = plt.subplots(figsize=(17.5, fig_h), facecolor=SURFACE)
        ax.set_facecolor(SURFACE)
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, fig_h)
        y = fig_h - 0.42

        s = seqs.get(orf, "")
        ax.text(0, y + 0.30, f"{orf} — 50 nearest proteins in ESM-2 650M embedding space",
                fontsize=14, color=INK, fontweight="bold", va="bottom")
        ax.text(0, y - 0.04,
                f"{len(s)} aa · ranked by cosine distance over 561,841 Swiss-Prot proteins · "
                f"DIAMOND e-value is for the same ORF↔neighbour pair",
                fontsize=9, color=INK_2, va="bottom")
        ax.text(0, y - 0.30, s if len(s) <= 118 else s[:117] + "…",
                fontsize=7.4, color=INK_MUTED, family="monospace", va="bottom")

        y -= 0.72
        x = 0.0
        for label, w in cols:
            ax.text(x, y, label, fontsize=8.6, color=INK_2, fontweight="bold", va="center")
            x += w
        y -= 0.10
        ax.plot([0, 1], [y, y], color=INK_2, lw=1.0)
        y -= 0.20

        for _, r in sub.iterrows():
            vals = [
                (str(int(r["rank"])), INK_MUTED, "normal"),
                (r.accession, INK_2, "normal"),
                (f"{r.cosine_distance:.4f}", INK, "normal"),
                ("no hit" if pd.isna(r.diamond_e) else f"{r.diamond_e:.1e}",
                 INK_MUTED if pd.isna(r.diamond_e) else INK, "normal"),
                (abbreviate(r.protein_name), INK, "normal"),
                (abbreviate(r.superfamily, 30) if isinstance(r.superfamily, str) else "—",
                 INK_2 if isinstance(r.superfamily, str) else INK_MUTED, "normal"),
                (abbreviate(r.organism, 26), INK_2, "normal"),
                (str(r.domain)[:4] if r.domain else "—", INK_2, "normal"),
                ("yes" if r.pathogen == "yes" else ("venom" if r.pathogen == "venom" else "—"),
                 HIT if r.pathogen == "yes" else (INK_2 if r.pathogen else INK_MUTED), "normal"),
                (abbreviate(r.vf, 13) if r.vf else "—", HIT if r.vf else INK_MUTED, "normal"),
            ]
            x = 0.0
            for (text, color, weight), (_, w) in zip(vals, cols):
                ax.text(x, y, text, fontsize=7.8, color=color, va="center", fontweight=weight)
                x += w
            y -= 0.213

        ax.plot([0, 1], [y + 0.10, y + 0.10], color=GRID, lw=0.8)
        nh = sub.diamond_e.isna().sum()
        ax.text(0, y - 0.12,
                f"DIAMOND (--very-sensitive, e≤1000): {len(sub) - nh}/{len(sub)} neighbours "
                f"align to {orf}; {nh} show no alignment at any threshold. "
                "High cosine similarity without alignment indicates compositional, not "
                "homologous, similarity.",
                fontsize=8, color=INK_2, va="top")
        ax.text(0, y - 0.40,
                "Path / VF describe the NEIGHBOUR protein and its source organism — they are "
                f"not attributes of {orf}. VF = DIAMOND homology to a VFDB curated virulence "
                "factor (e≤1e-5, ≥30% id).",
                fontsize=8, color=INK_MUTED, va="top")

        out = FIGDIR / f"{orf}_top50.png"
        fig.savefig(out, dpi=190, facecolor=SURFACE, bbox_inches="tight")
        plt.close(fig)
        print(f"  {out.name}: {len(sub)} rows, {nh} without alignment, "
              f"{(sub.pathogen == 'yes').sum()} pathogen-source, {(sub.vf != '').sum()} VF-homologous")


if __name__ == "__main__":
    d = build()
    print(f"annotated {len(d)} neighbour rows")
    render(d)
