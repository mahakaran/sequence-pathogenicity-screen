"""Characterise the AlphaFold 3 MSAs: how many DISTINCT proteins, and at what
identity/coverage. AF3 reports aligned segments, not independent relatives, so a
"depth" of 20 can be five proteins contributing repeated internal segments.

Usage:  python msa_depth.py <alphafold_extracted_dir> <out_csv>
"""
import glob
import math
import re
import sys
from collections import Counter
from pathlib import Path

SRC = Path(sys.argv[1])
OUT = Path(sys.argv[2])


def accession(header):
    """Bare accession, regardless of which spelling AF3 used.

    The unpaired file writes 'UniRef90_G4YP12/40-107'; the paired file writes
    'tr|G4YP12|G4YP12_PHYSP/40-107' for the same protein.
    """
    tok = header.lstrip(">").split()[0]
    tok = re.sub(r"/\d+-\d+.*$", "", tok)
    if "|" in tok:
        parts = tok.split("|")
        return parts[1] if len(parts) > 1 else tok
    return re.sub(r"^UniRef\d+_", "", tok)


def records(path):
    out = []
    for blk in open(path).read().split(">"):
        if blk.strip():
            h, *rest = blk.split("\n")
            out.append((h.strip(), "".join(rest)))
    return out


rows = []
for orf in ["orf1", "orf2", "orf3", "orf4", "orf5"]:
    chains = sorted({f.split("chains_")[1][0]
                     for f in glob.glob(str(SRC / orf / "msas" / "*.a3m"))})
    for chain in chains:
        files = sorted(glob.glob(str(SRC / orf / "msas" / f"*chains_{chain}.a3m")))
        srcs, idents, covs, segs, q = set(), [], [], 0, None
        for f in files:
            recs = records(f)
            if len(recs) < 2:
                continue
            q = recs[0][1]
            L = len(q)
            segs = max(segs, len(recs) - 1)
            for h, seq in recs[1:]:
                srcs.add(accession(h))
                core = [ch for ch in seq if not ch.islower()]
                pairs = [(a, b) for a, b in zip(q, core) if b != "-"]
                if not pairs:
                    continue
                idents.append(sum(a == b for a, b in pairs) / len(pairs) * 100)
                covs.append(len(pairs) / L * 100)
        if not srcs:
            rows.append({"orf": orf, "chain": chain, "segments": 0, "proteins": 0,
                         "top_source_share": None, "ident_min": None, "ident_max": None,
                         "cov_min": None, "cov_max": None, "query_entropy_bits": None})
            continue
        c = Counter(q)
        n = sum(c.values())
        H = -sum((v / n) * math.log2(v / n) for v in c.values())
        rows.append({
            "orf": orf, "chain": chain, "segments": segs, "proteins": len(srcs),
            "top_source_share": None,
            "ident_min": round(min(idents)), "ident_max": round(max(idents)),
            "cov_min": round(min(covs)), "cov_max": round(max(covs)),
            "query_entropy_bits": round(H, 2),
        })

import csv
FIELDS = ["orf", "chain", "segments", "proteins", "top_source_share",
          "ident_min", "ident_max", "cov_min", "cov_max", "query_entropy_bits"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r)
print(f"\nWrote {OUT}")
