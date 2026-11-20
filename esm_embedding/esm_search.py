"""Toy embedding-similarity search with ESM-2 (8M).

Embeds query proteins and a reference set with the same model, mean-pools
per-residue embeddings, and ranks references by cosine similarity.

References come from DS569k.parquet (Swiss-Prot, ~569k proteins). Two modes:

  proteinclip (default, fast)
      DS569k's `embedding` column = ESM-2 8M last layer, mean-pooled over
      residues, unit-normalised, then projected to 128-d by the pretrained
      ProteinCLIP adapter (github.com/wukevin/proteinclip). We apply the same
      pipeline to the queries and compare against the stored vectors directly,
      so no reference embedding is needed. (Verified: re-embedding Swiss-Prot
      entries this way reproduces the stored vectors with cosine = 1.0000.)

  esm
      Raw 320-d ESM-2 8M mean-pooled embeddings. References are re-embedded
      (~4 h for all of Swiss-Prot on an M2) and cached to disk.

Usage:
    python esm_search.py                          # ProteinCLIP vs all Swiss-Prot
    python esm_search.py --mode esm --limit 5000  # raw ESM-2 vs a random subset
    python esm_search.py --query my.fasta --top-k 25
"""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from Bio import SeqIO
from transformers import AutoModel, AutoTokenizer

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
MAX_RES = 1022  # ESM-2 positional limit (1024 tokens incl. <cls>/<eos>)
META_COLS = ["accession", "protein_name", "organism_name", "ncbi_taxonomy_phylum",
             "function", "sequence_length"]


def pick_device():
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


class Embedder:
    def __init__(self, model_name=MODEL_NAME, device=None):
        self.device = device or pick_device()
        self.tok = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name).eval().to(self.device)

    @torch.no_grad()
    def embed(self, seqs, max_tokens=16000, progress=False):
        """Mean-pooled, L2-normalised embeddings, in input order.

        Sequences are sorted by length and batched by token budget, which
        minimises padding and keeps memory bounded.
        """
        seqs = [s[:MAX_RES] for s in seqs]
        order = np.argsort([len(s) for s in seqs])
        dim = self.model.config.hidden_size
        out = np.zeros((len(seqs), dim), dtype=np.float32)

        batches, cur, cur_max = [], [], 0
        for i in order:
            L = len(seqs[i]) + 2
            if cur and max(cur_max, L) * (len(cur) + 1) > max_tokens:
                batches.append(cur)
                cur, cur_max = [], 0
            cur.append(i)
            cur_max = max(cur_max, L)
        if cur:
            batches.append(cur)

        t0, done = time.time(), 0
        for b, idx in enumerate(batches):
            enc = self.tok([seqs[i] for i in idx], return_tensors="pt",
                           padding=True).to(self.device)
            h = self.model(**enc).last_hidden_state            # (B, T, D)
            mask = enc["attention_mask"].clone()
            lengths = mask.sum(1)                              # residues + 2
            mask[:, 0] = 0                                     # drop <cls>
            mask[torch.arange(len(idx)), lengths - 1] = 0      # drop <eos>
            mask = mask.unsqueeze(-1).to(h.dtype)
            pooled = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
            pooled = torch.nn.functional.normalize(pooled, dim=-1)
            out[idx] = pooled.float().cpu().numpy()

            done += len(idx)
            if progress and (b % 50 == 0 or done == len(seqs)):
                rate = done / max(time.time() - t0, 1e-9)
                eta = (len(seqs) - done) / max(rate, 1e-9)
                print(f"  {done:>7}/{len(seqs)}  {rate:6.0f} seq/s  ETA {eta/60:5.1f} min",
                      flush=True)
        return out


class ProteinCLIP:
    """Pretrained ProteinCLIP adapter for ESM-2 8M: 320-d -> 128-d unit vector."""

    def __init__(self, path=HERE / "pretrained" / "proteinclip_esm2_6.onnx"):
        import onnxruntime as ort
        self.sess = ort.InferenceSession(str(path), providers=["CPUExecutionProvider"])

    def __call__(self, x):
        x = x / np.linalg.norm(x, axis=1, keepdims=True)
        y = self.sess.run(None, {"input": x.astype(np.float32)})[0]
        return y / np.linalg.norm(y, axis=1, keepdims=True)


def load_references(path, limit=None, seed=0, with_embedding=False):
    cols = META_COLS + ["sequence"] + (["embedding"] if with_embedding else [])
    meta = pq.read_table(path, columns=cols).to_pandas()
    if limit and limit < len(meta):
        meta = meta.sample(n=limit, random_state=seed).reset_index(drop=True)
    return meta


def reference_embeddings(embedder, meta, cache_dir, tag):
    cache = cache_dir / f"ref_{tag}.npy"
    ids = cache_dir / f"ref_{tag}_accessions.txt"
    if cache.exists() and ids.exists():
        cached_ids = ids.read_text().split()
        if cached_ids == meta["accession"].tolist():
            print(f"Loaded cached reference embeddings: {cache.name}")
            return np.load(cache)
        print("Cache does not match reference set; re-embedding.")
    print(f"Embedding {len(meta)} reference sequences on {embedder.device}...")
    emb = embedder.embed(meta["sequence"].tolist(), progress=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.save(cache, emb)
    ids.write_text("\n".join(meta["accession"]))
    return emb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", type=Path, default=PROJECT / "screening-exercise-orfs.fasta",
                    help="protein FASTA of query sequences")
    ap.add_argument("--refs", type=Path, default=PROJECT / "DS569k.parquet",
                    help="parquet with reference sequences + annotations")
    ap.add_argument("--mode", choices=["proteinclip", "esm"], default="proteinclip",
                    help="proteinclip: use DS569k's stored 128-d vectors; "
                         "esm: raw 320-d ESM-2, re-embed references")
    ap.add_argument("--limit", type=int, default=None,
                    help="random subset of N references (for quick tests)")
    ap.add_argument("--top-k", type=int, default=10)
    ap.add_argument("--min-len", type=int, default=0,
                    help="skip query ORFs shorter than this")
    ap.add_argument("--out-dir", type=Path, default=HERE / "results")
    ap.add_argument("--cache-dir", type=Path, default=HERE / "cache")
    args = ap.parse_args()

    queries = [r for r in SeqIO.parse(args.query, "fasta")
               if len(r.seq) >= args.min_len]
    q_seqs = [str(r.seq).rstrip("*").upper() for r in queries]
    print(f"{len(queries)} query sequences from {args.query.name}")

    embedder = Embedder()
    print(f"Model: {MODEL_NAME}  device: {embedder.device}")

    subset = "all" if not args.limit else f"n{args.limit}"
    if args.mode == "proteinclip":
        meta = load_references(args.refs, args.limit, with_embedding=True)
        ref_emb = np.vstack(meta.pop("embedding").values).astype(np.float32)
        ref_emb /= np.linalg.norm(ref_emb, axis=1, keepdims=True)
        q_emb = ProteinCLIP()(embedder.embed(q_seqs))
        tag = f"proteinclip_esm2_t6_8M_{subset}"
        print(f"Using {len(meta)} stored ProteinCLIP reference embeddings")
    else:
        meta = load_references(args.refs, args.limit)
        tag = f"esm2_t6_8M_{subset}"
        ref_emb = reference_embeddings(embedder, meta, args.cache_dir, tag)
        q_emb = embedder.embed(q_seqs)

    sims = q_emb @ ref_emb.T                                   # cosine similarity
    rows = []
    for qi, rec in enumerate(queries):
        s = sims[qi]
        mu, sd = s.mean(), s.std()
        top = np.argsort(-s)[: args.top_k]
        print(f"\n=== {rec.id}  (len {len(q_seqs[qi])})  "
              f"background cos mean={mu:.3f} sd={sd:.3f}")
        print(f"{'rank':>4} {'cos':>6} {'z':>6}  {'accession':<10} {'len':>5}  protein | organism")
        for rank, j in enumerate(top, 1):
            m = meta.iloc[j]
            z = (s[j] - mu) / sd
            print(f"{rank:>4} {s[j]:6.3f} {z:6.2f}  {m.accession:<10} {m.sequence_length:>5}  "
                  f"{str(m.protein_name)[:50]} | {str(m.organism_name)[:40]}")
            rows.append({"query": rec.id, "query_len": len(q_seqs[qi]), "rank": rank,
                         "cosine": float(s[j]), "z_score": float(z),
                         **{c: m[c] for c in META_COLS}})

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / f"top{args.top_k}_{tag}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    np.save(args.out_dir / "query_embeddings.npy", q_emb)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
