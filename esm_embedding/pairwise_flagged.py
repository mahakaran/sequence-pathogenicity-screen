"""Rigorous pairwise alignment: each ORF vs its flagged neighbours.

Flagged = source organism is a pathogen, or venomous/toxin-producing, or the
neighbour has DIAMOND homology to a VFDB curated virulence factor.

DIAMOND is a fast heuristic and can miss weak similarity in short, low-complexity
sequences. This runs full Smith-Waterman (local, BLOSUM62) on every flagged pair
and calibrates significance empirically: for each pair, the ORF is shuffled N times
and realigned, giving a null distribution of scores. The empirical p-value is the
fraction of shuffles scoring at least as high as the real alignment.

Shuffling preserves amino-acid composition exactly, so a pair that scores highly
only because both sequences are, say, hydrophobic will also score highly against
the shuffles and yield a non-significant p-value. That is precisely the
compositional-vs-homologous distinction we need.

Usage:  python pairwise_flagged.py [--shuffles 2000]
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
OUT = RES / "flagged_alignments.parquet"


def aligner():
    a = Align.PairwiseAligner(scoring="blastp")   # BLOSUM62, gap -11/-1
    a.mode = "local"
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shuffles", type=int, default=2000)
    args = ap.parse_args()

    d = pd.read_parquet(RES / "top50_annotated.parquet")
    flagged = d[(d.pathogen != "") | (d.vf != "")].copy()

    meta = pq.read_table(PROJECT / "DS569k.parquet",
                         columns=["accession", "sequence"]).to_pandas().set_index("accession")
    seqs = {}
    for blk in (PROJECT / "screening-exercise-orfs.fasta").read_text().split(">"):
        if blk.strip():
            h, *rest = blk.strip().split("\n")
            seqs[h.split("|")[-1].split(":")[0]] = "".join(rest).rstrip("*").upper()

    al = aligner()
    rng = np.random.default_rng(0)
    rows = []
    for i, (_, r) in enumerate(flagged.iterrows(), 1):
        q = seqs[r.orf]
        s = str(meta.loc[r.accession, "sequence"]).upper()
        # Biopython rejects non-standard residues against BLOSUM62
        s = "".join(c for c in s if c in "ACDEFGHIKLMNPQRSTVWY")
        q2 = "".join(c for c in q if c in "ACDEFGHIKLMNPQRSTVWY")
        if not s or not q2:
            continue
        real = al.score(q2, s)

        arr = np.array(list(q2))
        null = np.empty(args.shuffles)
        for k in range(args.shuffles):
            rng.shuffle(arr)
            null[k] = al.score("".join(arr), s)
        # +1 smoothing: with N shuffles the smallest reportable p-value is 1/(N+1)
        p = (np.sum(null >= real) + 1) / (args.shuffles + 1)

        aln = al.align(q2, s)[0]
        a_q, a_s = str(aln).split("\n")[0], str(aln).split("\n")[2]
        matches = sum(1 for x, y in zip(a_q, a_s) if x == y and x != "-")
        alen = sum(1 for x, y in zip(a_q, a_s) if x != "-" and y != "-")

        rows.append({
            "orf": r.orf, "rank": int(r["rank"]), "accession": r.accession,
            "protein_name": r.protein_name, "organism": r.organism, "domain": r.domain,
            "pathogen": r.pathogen, "vf": r.vf,
            "cosine_distance": r.cosine_distance,
            "sw_score": float(real),
            "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "z": float((real - null.mean()) / null.std()) if null.std() else np.nan,
            "p_empirical": float(p),
            "pct_identity": 100.0 * matches / alen if alen else np.nan,
            "aln_length": int(alen),
        })
        if i % 20 == 0:
            print(f"  {i}/{len(flagged)} pairs aligned", flush=True)

    out = pd.DataFrame(rows)
    out.to_parquet(OUT, index=False)

    print(f"\n{len(out)} flagged pairs aligned, {args.shuffles} shuffles each")
    print(f"min empirical p-value observed: {out.p_empirical.min():.4f} "
          f"(floor is {1 / (args.shuffles + 1):.5f})")
    sig = out[out.p_empirical < 0.05]
    print(f"pairs with p < 0.05: {len(sig)}")
    if len(sig):
        print("\n  " + "-" * 94)
        for _, r in sig.sort_values("p_empirical").iterrows():
            print(f"  {r.orf} #{r['rank']:<3} {r.accession}  SW={r.sw_score:.0f}  "
                  f"z={r.z:+.2f}  p={r.p_empirical:.4f}  id={r.pct_identity:.0f}%  "
                  f"{str(r.protein_name)[:38]}")
    print("\nPer-ORF summary (z-score of real alignment vs shuffled null):")
    for o, g in out.groupby("orf"):
        print(f"  {o}: n={len(g):<3} max z={g.z.max():+.2f}  min p={g.p_empirical.min():.4f}  "
              f"max identity={g.pct_identity.max():.0f}%")


if __name__ == "__main__":
    main()
