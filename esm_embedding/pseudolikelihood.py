"""Are the ORFs protein-like? ESM-2 pseudo-log-likelihood against matched controls.

PLL = (1/L) * sum_i log P(x_i | x_without_i), masking one position at a time.
Reported per residue so sequences of different length stay comparable; even so,
PLL still varies systematically with length, which is why every control here is
length-matched to the ORF it is compared against.

Controls, per ORF (n each):
  uniprot     real Swiss-Prot proteins of the same length (+-10%)        [positive]
  shuffled    those same proteins, residues shuffled                     [negative]
  comp_match  the ORF's OWN residues shuffled - same composition,        [negative]
              different order. The strictest control: it asks whether
              the ORDER is protein-like, not whether the composition is.
  uniform     uniform random over the 20 amino acids                     [negative]

Usage:
    python pseudolikelihood.py --bench            # time a few sequences first
    python pseudolikelihood.py --n 100
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
RES = HERE / "results"
MODEL_DIR = HERE / "model_650m"
OUT = RES / "pll_scores.parquet"
AA = "ACDEFGHIKLMNPQRSTVWY"
GROUPS = ["uniprot", "shuffled", "comp_match", "uniform"]


def load():
    tok = AutoTokenizer.from_pretrained(str(MODEL_DIR))
    model = AutoModelForMaskedLM.from_pretrained(str(MODEL_DIR)).eval()
    dev = "mps" if torch.backends.mps.is_available() else "cpu"
    return tok, model.to(dev), dev


@torch.no_grad()
def pll(seq, tok, model, dev, batch=32):
    """Mean per-residue pseudo-log-likelihood."""
    enc = tok([seq], return_tensors="pt")
    ids = enc["input_ids"][0].to(dev)
    L = len(seq)
    positions = torch.arange(1, L + 1, device=dev)       # skip <cls> and <eos>
    total = 0.0
    for s in range(0, L, batch):
        idx = positions[s:s + batch]
        rep = ids.unsqueeze(0).repeat(len(idx), 1)
        rep[torch.arange(len(idx)), idx] = tok.mask_token_id
        out = model(input_ids=rep, attention_mask=torch.ones_like(rep)).logits
        lp = torch.log_softmax(out[torch.arange(len(idx)), idx], dim=-1)
        total += lp[torch.arange(len(idx)), ids[idx]].sum().item()
    return total / L


def build_controls(orfs, n, rng):
    """Length-matched control sequences for every ORF."""
    t = pq.read_table(PROJECT / "DS569k.parquet",
                      columns=["accession", "sequence"]).to_pandas()
    t["L"] = t.sequence.str.len()
    clean = t[t.sequence.str.fullmatch(f"[{AA}]+")]
    rows = []
    for name, seq in orfs.items():
        L = len(seq)
        lo, hi = int(L * 0.9), int(np.ceil(L * 1.1))
        pool = clean[(clean.L >= lo) & (clean.L <= hi)]
        if len(pool) < n:                      # widen if the band is sparse
            lo, hi = int(L * 0.75), int(np.ceil(L * 1.35))
            pool = clean[(clean.L >= lo) & (clean.L <= hi)]
        take = pool.sample(min(n, len(pool)), random_state=int(rng.integers(1e9)))
        for acc, s in zip(take.accession, take.sequence):
            rows.append({"orf": name, "group": "uniprot", "label": acc,
                         "sequence": s, "length": len(s)})
            sh = list(s)
            rng.shuffle(sh)
            rows.append({"orf": name, "group": "shuffled", "label": f"shuf_{acc}",
                         "sequence": "".join(sh), "length": len(s)})
        for i in range(n):
            c = list(seq)
            rng.shuffle(c)
            rows.append({"orf": name, "group": "comp_match", "label": f"comp{i}",
                         "sequence": "".join(c), "length": L})
            rows.append({"orf": name, "group": "uniform", "label": f"unif{i}",
                         "sequence": "".join(rng.choice(list(AA), size=L)), "length": L})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--bench", action="store_true")
    args = ap.parse_args()

    orfs = {}
    for blk in (PROJECT / "screening-exercise-orfs.fasta").read_text().split(">"):
        if blk.strip():
            h, *rest = blk.strip().split("\n")
            s = "".join(rest).rstrip("*").upper()
            orfs[h.split("|")[-1].split(":")[0]] = "".join(c for c in s if c in AA)

    tok, model, dev = load()
    print(f"device {dev}")

    if args.bench:
        import time
        for name, s in list(orfs.items())[:2]:
            t0 = time.time()
            v = pll(s, tok, model, dev)
            print(f"  {name} ({len(s)} aa): PLL {v:+.4f}  [{time.time() - t0:.1f}s]")
        return

    rng = np.random.default_rng(0)
    ctl = build_controls(orfs, args.n, rng)
    print(f"controls: {len(ctl):,} sequences "
          f"({args.n} per group x {len(GROUPS)} groups x {len(orfs)} ORFs)")
    print(f"median control length by ORF: "
          f"{ctl.groupby('orf').length.median().to_dict()}")

    rows = []
    for name, s in orfs.items():
        rows.append({"orf": name, "group": "ORF", "label": name,
                     "sequence": s, "length": len(s), "pll": pll(s, tok, model, dev)})
        print(f"  {name}: PLL {rows[-1]['pll']:+.4f}", flush=True)

    for i, (_, r) in enumerate(ctl.iterrows(), 1):
        rows.append({**r.to_dict(), "pll": pll(r.sequence, tok, model, dev)})
        if i % 100 == 0:
            print(f"  controls {i:,}/{len(ctl):,}", flush=True)

    d = pd.DataFrame(rows).drop(columns=["sequence"])
    d.to_parquet(OUT, index=False)
    print(f"\nWrote {OUT}")
    print(d.groupby("group").pll.describe()[["count", "mean", "std", "min", "max"]])


if __name__ == "__main__":
    main()
