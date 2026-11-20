"""Embed the Swiss-Prot proteins absent from embeddings_esm2_650m.h5.

Uses the pooling convention verified against that file (mean over residue
positions, <cls>/<eos> excluded), so the new vectors are directly comparable.

Writes results/missing_650m.npy + results/missing_650m_accessions.txt, saving
partial progress so an interrupted run resumes instead of restarting.

Usage:
    python embed_missing.py --limit 200      # benchmark
    python embed_missing.py                  # full run
"""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
SRC = RES / "missing_to_embed.parquet"
OUT = RES / "missing_650m.npy"
ACC = RES / "missing_650m_accessions.txt"
CKPT = RES / "missing_650m_partial.npz"          # legacy, migrated on resume
CKPT_BIN = RES / "missing_650m_partial.f32"
CKPT_ACC = RES / "missing_650m_partial_acc.txt"
MODEL_DIR = HERE / "model_650m"
MAX_RES = 1022
DIM = 1280


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--max-tokens", type=int, default=12000,
                    help="token budget per batch (length x batch size)")
    ap.add_argument("--checkpoint-every", type=int, default=2000)
    args = ap.parse_args()

    df = pd.read_parquet(SRC)
    if args.limit:
        df = df.head(args.limit)
    seqs = [s[:MAX_RES].upper() for s in df.sequence]
    accs = list(df.accession)

    # Checkpoint is append-only: one float32 vector per row in a raw .f32 file with a
    # parallel accession list. Rewriting a single npz each time (vstack over everything
    # embedded so far) grows quadratically and thrashes memory on a long run.
    # Only the SET of finished accessions is kept in memory; the vectors themselves stay
    # on disk and are read back once at the end. Holding 28k x 1280 floats in a dict
    # alongside the 2.5 GB model exhausts physical memory and the run swaps itself dead.
    done = {}
    if CKPT_BIN.exists() and CKPT_ACC.exists() and not args.limit:
        prev = CKPT_ACC.read_text().split()
        n_rows = CKPT_BIN.stat().st_size // (DIM * 4)
        n = min(len(prev), n_rows)             # tolerate a half-written tail
        done = {a: None for a in prev[:n]}
        print(f"resuming: {n:,} already embedded (vectors stay on disk)")
    elif CKPT.exists() and not args.limit:
        z = np.load(CKPT, allow_pickle=True)   # migrate the old npz checkpoint
        done = dict(zip(z["acc"].tolist(), z["emb"]))
        with open(CKPT_BIN, "wb") as fh:
            np.vstack([done[a] for a in done]).astype(np.float32).tofile(fh)
        CKPT_ACC.write_text("\n".join(done))
        print(f"migrated {len(done):,} embeddings from the old checkpoint format")
    todo = [i for i, a in enumerate(accs) if a not in done]

    tok = AutoTokenizer.from_pretrained(str(MODEL_DIR))
    model = AutoModel.from_pretrained(str(MODEL_DIR)).eval()
    dev = "mps" if torch.backends.mps.is_available() else "cpu"
    model.to(dev)
    print(f"device {dev}; {len(todo):,} sequences to embed")

    # batch by token budget, shortest first, to minimise padding
    todo.sort(key=lambda i: len(seqs[i]))
    batches, cur, cur_max = [], [], 0
    for i in todo:
        L = len(seqs[i]) + 2
        if cur and max(cur_max, L) * (len(cur) + 1) > args.max_tokens:
            batches.append(cur)
            cur, cur_max = [], 0
        cur.append(i)
        cur_max = max(cur_max, L)
    if cur:
        batches.append(cur)

    t0, n = time.time(), 0
    ck_bin = open(CKPT_BIN, "ab") if not args.limit else None
    ck_acc = open(CKPT_ACC, "a") if not args.limit else None
    with torch.no_grad():
        for b, idx in enumerate(batches):
            enc = tok([seqs[i] for i in idx], return_tensors="pt", padding=True).to(dev)
            h = model(**enc).last_hidden_state
            mask = enc["attention_mask"].clone()
            lengths = mask.sum(1)
            mask[:, 0] = 0                                   # <cls>
            mask[torch.arange(len(idx)), lengths - 1] = 0    # <eos>
            mask = mask.unsqueeze(-1).to(h.dtype)
            pooled = ((h * mask).sum(1) / mask.sum(1).clamp(min=1)).float().cpu().numpy()
            for i in idx:
                done[accs[i]] = None                     # track completion, not the vector
            if ck_bin is not None:                       # append this batch only
                pooled.astype(np.float32).tofile(ck_bin)
                ck_acc.write("".join(f"{accs[i]}\n" for i in idx))
                ck_bin.flush()
                ck_acc.flush()
            n += len(idx)
            if b % 20 == 0 or n == len(todo):
                rate = n / max(time.time() - t0, 1e-9)
                print(f"  {n:>6,}/{len(todo):,}  {rate:5.1f} seq/s  "
                      f"ETA {(len(todo) - n) / max(rate, 1e-9) / 60:5.1f} min", flush=True)
    if ck_bin is not None:
        ck_bin.close()
        ck_acc.close()

    # Read the checkpoint back and reorder it to match the source parquet.
    order = CKPT_ACC.read_text().split()
    raw = np.fromfile(CKPT_BIN, dtype=np.float32).reshape(-1, DIM)
    pos = {a: i for i, a in enumerate(order[:len(raw)])}
    keep = [a for a in accs if a in pos]
    emb = raw[[pos[a] for a in keep]]
    np.save(OUT, emb)
    ACC.write_text("\n".join(keep))
    print(f"\nWrote {OUT} {emb.shape} and {ACC.name}")
    if len(keep) != len(accs):
        print(f"NOTE: {len(accs) - len(keep):,} sequences were not embedded")
    if not args.limit:
        for p in (CKPT, CKPT_BIN, CKPT_ACC):
            p.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
