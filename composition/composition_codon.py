"""GC/AT composition and codon-usage analysis of the screening sequence.

Signals computed:
  1. Base composition: GC, skews, sliding-window GC, GC1/GC2/GC3 per open frame
  2. Dinucleotide observed/expected (CpG, TpA, ...)
  3. Codon usage of the two full-length open frames (+1 = ORF1, -1 = ORF4) vs a
     reference panel (Kazusa codon usage DB): CAI, synonymous-codon
     log-likelihood, GC3 distance, ENC
  4. "Designed?" checks: restriction sites, repeats, homopolymers, and how
     surprising it is to have long open frames on BOTH strands (shuffle tests)

Usage:  python composition_codon.py      (writes report.txt next to this file)
"""

import itertools
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
SEQ_FILE = HERE.parent / "sequence.fasta"
REF_DIR = HERE / "ref_tables"
RNG = np.random.default_rng(0)

# Reference panel: (Kazusa taxid, label, group, approx genome GC %)
PANEL = [
    (83333, "Escherichia coli K-12", "Bacteria", 51),
    (1423, "Bacillus subtilis", "Bacteria", 43),
    (158879, "Staphylococcus aureus N315", "Bacteria", 33),
    (195102, "Clostridium perfringens 13", "Bacteria", 28),
    (1902, "Streptomyces coelicolor", "Bacteria", 72),
    (9606, "Homo sapiens", "Eukaryota", 41),
    (4932, "Saccharomyces cerevisiae", "Eukaryota", 38),
    (5476, "Candida albicans", "Eukaryota", 33),
    (3702, "Arabidopsis thaliana", "Eukaryota", 36),
    (7227, "Drosophila melanogaster", "Eukaryota", 42),
    (6239, "Caenorhabditis elegans", "Eukaryota", 35),
    (44689, "Dictyostelium discoideum", "Eukaryota", 22),
    (36329, "Plasmodium falciparum 3D7", "Eukaryota", 19),
    (10254, "Vaccinia virus WR", "Virus", 33),
    (10665, "Enterobacteria phage T4", "Virus", 35),
    (46015, "Autographa californica NPV", "Virus", 41),
    (10298, "Human herpesvirus 1", "Virus", 68),
]

STD = CodonTable.unambiguous_dna_by_id[1]
CODONS = ["".join(c) for c in itertools.product("TCAG", repeat=3)]
AA_OF = {c: STD.forward_table.get(c, "*") for c in CODONS}
SYN = {}
for c, a in AA_OF.items():
    SYN.setdefault(a, []).append(c)

RESTRICTION_SITES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "XhoI": "CTCGAG",
    "NdeI": "CATATG", "NcoI": "CCATGG", "XbaI": "TCTAGA", "BglII": "AGATCT",
    "NotI": "GCGGCCGC", "SalI": "GTCGAC", "KpnI": "GGTACC", "SacI": "GAGCTC",
    "PstI": "CTGCAG", "SpeI": "ACTAGT", "NheI": "GCTAGC", "EcoRV": "GATATC",
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "SapI": "GCTCTTC", "AscI": "GGCGCGCC",
}


# ---------------------------------------------------------------- references
def fetch_table(taxid):
    """Kazusa codon table -> {codon: count}. Cached in ref_tables/."""
    REF_DIR.mkdir(exist_ok=True)
    f = REF_DIR / f"{taxid}.txt"
    if not f.exists():
        url = f"http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={taxid}&aa=1&style=N"
        html = urllib.request.urlopen(url, timeout=60).read().decode("latin-1")
        f.write_text(html)
    html = f.read_text()
    pre = re.search(r"<PRE>(.*?)</PRE>", html, re.S | re.I)
    if not pre:
        raise ValueError(f"no codon table for taxid {taxid}")
    counts = {}
    # NB: Kazusa columns can run together ("AAU N 0.86123.5 (369612)"), so match
    # loosely up to the count in parentheses rather than on whitespace.
    for cod, cnt in re.findall(r"([ACGU]{3}) [A-Z*] [^()]*\(\s*(\d+)\)", pre.group(1)):
        counts[cod.replace("U", "T")] = int(cnt)
    if len(counts) != 64:
        raise ValueError(f"parsed {len(counts)} codons for taxid {taxid}")
    n_codons = re.search(r"([\d,]+)\s+codons", html)
    return counts, (int(n_codons.group(1).replace(",", "")) if n_codons else sum(counts.values()))


def syn_freqs(counts, pseudo=0.5):
    """P(codon | amino acid) with a small pseudocount."""
    p = {}
    for aa, cods in SYN.items():
        tot = sum(counts[c] + pseudo for c in cods)
        for c in cods:
            p[c] = (counts[c] + pseudo) / tot
    return p


def gc3_of_table(counts):
    sense = [c for c in CODONS if AA_OF[c] not in "*MW"]
    tot = sum(counts[c] for c in sense)
    return 100 * sum(counts[c] for c in sense if c[2] in "GC") / tot


# ---------------------------------------------------------------- metrics
def gc(s):
    return 100 * sum(b in "GC" for b in s) / len(s)


def codons_of(s):
    return [str(s[i:i + 3]) for i in range(0, len(s) - len(s) % 3, 3)]


def cai(codons, counts):
    w = {}
    for aa, cods in SYN.items():
        m = max(counts[c] for c in cods)
        for c in cods:
            w[c] = max(counts[c], 0.5) / max(m, 0.5)
    use = [c for c in codons if AA_OF[c] not in "*MW"]
    return math.exp(sum(math.log(w[c]) for c in use) / len(use))


def syn_loglik(codons, counts):
    """Mean log2 P(codon | aa, organism): which organism's codon *choice* fits best."""
    p = syn_freqs(counts)
    use = [c for c in codons if len(SYN[AA_OF[c]]) > 1 and AA_OF[c] != "*"]
    return sum(math.log2(p[c]) for c in use) / len(use)


def enc_wright(codons):
    """Wright's effective number of codons (20 = max bias, 61 = uniform)."""
    by_aa = Counter(codons)
    F = {2: [], 3: [], 4: [], 6: []}
    for aa, cods in SYN.items():
        k = len(cods)
        if aa == "*" or k == 1:
            continue
        n = sum(by_aa[c] for c in cods)
        if n < 2:
            continue
        f = (n * sum((by_aa[c] / n) ** 2 for c in cods) - 1) / (n - 1)
        F[k].append(f)
    enc = 2.0
    for k, nk in [(2, 9), (3, 1), (4, 5), (6, 3)]:
        vals = [v for v in F[k] if v > 0]
        if vals:
            enc += nk / np.mean(vals)
        elif k == 3:
            enc += 1 / np.mean([np.mean(F[2] or [0.5]), np.mean(F[4] or [0.25])])
    return min(enc, 61.0)


def dinuc_oe(s):
    n = len(s)
    mono = Counter(s)
    di = Counter(s[i:i + 2] for i in range(n - 1))
    out = {}
    for a, b in itertools.product("ACGT", repeat=2):
        exp = mono[a] * mono[b] / n
        out[a + b] = di[a + b] / exp if exp else float("nan")
    return out


def open_frames(s):
    """Longest stop-free run (codons) in each of the 6 frames."""
    res = {}
    for strand, t in (("+", s), ("-", str(Seq(s).reverse_complement()))):
        for f in range(3):
            cods = codons_of(t[f:])
            run = best = 0
            for c in cods:
                run = 0 if AA_OF[c] == "*" else run + 1
                best = max(best, run)
            res[f"{strand}{f + 1}"] = (best, len(cods))
    return res


def fully_open(s):
    return {k: best == n for k, (best, n) in open_frames(s).items()}


# ---------------------------------------------------------------- main
def main():
    out = []
    say = out.append
    nt = str(next(SeqIO.parse(SEQ_FILE, "fasta")).seq).upper()
    L = len(nt)
    frames = {"+1 (ORF1)": nt[: L // 3 * 3],
              "-1 (ORF4, SOD-like)": str(Seq(nt).reverse_complement())[: L // 3 * 3]}

    # 1. composition ---------------------------------------------------------
    comp = Counter(nt)
    say(f"Sequence length: {L} nt")
    say("Base counts: " + "  ".join(f"{b}={comp[b]} ({100 * comp[b] / L:.1f}%)" for b in "ACGT"))
    say(f"GC = {gc(nt):.1f}%   AT = {100 - gc(nt):.1f}%")
    say(f"GC skew (G-C)/(G+C) = {(comp['G'] - comp['C']) / (comp['G'] + comp['C']):+.2f}   "
        f"AT skew (A-T)/(A+T) = {(comp['A'] - comp['T']) / (comp['A'] + comp['T']):+.2f}")
    win = [gc(nt[i:i + 60]) for i in range(0, L - 59, 15)]
    say(f"Sliding GC (60 nt windows, step 15): min {min(win):.0f}%  max {max(win):.0f}%  -> "
        + " ".join(f"{w:.0f}" for w in win))

    say("\nGC by codon position (GC1 / GC2 / GC3):")
    for name, s in frames.items():
        cods = codons_of(s)
        g = [100 * sum(c[i] in "GC" for c in cods) / len(cods) for i in range(3)]
        say(f"  {name:<22} {g[0]:5.1f} / {g[1]:5.1f} / {g[2]:5.1f}   ({len(cods)} codons)")

    # 2. dinucleotides -------------------------------------------------------
    oe = dinuc_oe(nt)
    say("\nDinucleotide observed/expected (1.0 = as expected from base composition):")
    say("  " + "  ".join(f"{k}={v:.2f}" for k, v in oe.items()))
    say(f"  CpG o/e = {oe['CG']:.2f}  (vertebrates & many vertebrate viruses ~0.2-0.4; "
        "most bacteria/phage ~0.8-1.2)")
    say(f"  TpA o/e = {oe['TA']:.2f}  (typically depleted, ~0.6-0.8, in most genomes)")

    # 3. codon usage vs panel -----------------------------------------------
    tables = {}
    for taxid, label, group, ggc in PANEL:
        try:
            tables[taxid] = fetch_table(taxid)
        except Exception as e:  # noqa: BLE001
            say(f"  (could not load {label}: {e})")

    for name, s in frames.items():
        cods = codons_of(s)
        g3 = 100 * sum(c[2] in "GC" for c in cods if AA_OF[c] not in "*MW") / \
            sum(AA_OF[c] not in "*MW" for c in cods)
        say(f"\nCodon usage, frame {name}: ENC = {enc_wright(cods):.1f} "
            f"(20 = extreme bias, 61 = none; noisy at {len(cods)} codons), GC3 = {g3:.1f}%")
        rows = []
        for taxid, label, group, ggc in PANEL:
            if taxid not in tables:
                continue
            counts, n = tables[taxid]
            rows.append((syn_loglik(cods, counts), cai(cods, counts), gc3_of_table(counts),
                         label, group, n))
        rows.sort(reverse=True)
        say(f"  {'rank':>4}  {'synLL':>6}  {'CAI':>5}  {'refGC3':>6}  organism (group, codons in table)")
        for i, (ll, c, rg3, label, group, n) in enumerate(rows, 1):
            say(f"  {i:>4}  {ll:6.2f}  {c:5.2f}  {rg3:5.1f}%  {label} ({group}, {n:,})")
        say("  synLL = mean log2 P(codon | amino acid) under that organism's usage "
            "(higher = codon choices fit better). CAI = codon adaptation index.")

    # 4. designed? -----------------------------------------------------------
    say("\nOpen reading frames (longest stop-free run / frame length, codons):")
    of = open_frames(nt)
    say("  " + "  ".join(f"{k}: {b}/{n}" for k, (b, n) in of.items()))
    stops = Counter(c for s in frames.values() for c in codons_of(s) if AA_OF[c] == "*")
    say(f"  Stops in +1/-1 frames: {dict(stops) or 'none'}")

    # (a) mononucleotide shuffles
    N = 20000
    arr = np.array(list(nt))
    both = any_opp = 0
    for _ in range(N):
        RNG.shuffle(arr)
        fo = fully_open("".join(arr))
        both += fo["+1"] and fo["-1"]
        any_opp += any(fo[p] for p in ("+1", "+2", "+3")) and any(fo[m] for m in ("-1", "-2", "-3"))
    say(f"\nShuffle test A - random sequences with the same base composition (n={N:,}):")
    say(f"  P(+1 and -1 both stop-free over all 96 codons) = {both / N:.5f}")
    say(f"  P(any frame on each strand fully stop-free)     = {any_opp / N:.5f}")

    # (b) hold the ORF4 protein fixed, re-draw synonymous codons from each reference
    #     organism's usage, and ask how often the opposite strand stays open.
    orf4 = codons_of(frames["-1 (ORF4, SOD-like)"])
    rc_is_stop = {c for c in CODONS
                  if AA_OF[str(Seq(c).reverse_complement())] == "*"}
    used = Counter(orf4)
    say(f"\nCodons whose reverse complement is a stop codon ({', '.join(sorted(rc_is_stop))}): "
        f"used {sum(used[c] for c in rc_is_stop)} times in ORF4's 96 codons")
    for aa in "LSF":
        say(f"  {aa}: " + ", ".join(f"{c}={used[c]}" + ("*" if c in rc_is_stop else "")
                                    for c in sorted(SYN[aa])))

    say(f"\nShuffle test B - ORF4's protein kept, synonymous codons re-drawn from each "
        f"organism's usage (n={N // 4:,} each):")
    say(f"  {'P(opposite strand stop-free)':>30}  organism")
    for taxid, label, group, ggc in PANEL:
        if taxid not in tables:
            continue
        pr = syn_freqs(tables[taxid][0])
        pick = {aa: (cods, np.array([pr[c] for c in cods]) /
                     sum(pr[c] for c in cods)) for aa, cods in SYN.items()}
        hits = 0
        for _ in range(N // 4):
            new = [pick[AA_OF[c]][0][RNG.choice(len(pick[AA_OF[c]][0]), p=pick[AA_OF[c]][1])]
                   for c in orf4]
            anti = codons_of(str(Seq("".join(new)).reverse_complement()))
            hits += all(AA_OF[c] != "*" for c in anti)
        say(f"  {hits / (N // 4):>30.5f}  {label} ({group})")

    say("\nRestriction sites (both strands, palindromic sites count once):")
    found = []
    rc = str(Seq(nt).reverse_complement())
    for enz, site in RESTRICTION_SITES.items():
        pos = [m.start() + 1 for m in re.finditer(f"(?={site})", nt)]
        rsite = str(Seq(site).reverse_complement())
        if rsite != site:
            pos += [L - m.start() - len(site) + 1 for m in re.finditer(f"(?={site})", rc)]
        if pos:
            found.append(f"{enz} ({site}) at {sorted(pos)}")
    say("  " + ("; ".join(found) if found else "none of the common cloning sites"))
    exp6 = (L - 5) * np.prod([comp[b] / L for b in "AGATCT"])
    say(f"  (for reference: expected count of a given AT-rich 6-mer like AGATCT here ~{exp6:.2f})")

    runs = [(m.group(0)[0], len(m.group(0)), m.start() + 1) for m in re.finditer(r"(A+|C+|G+|T+)", nt)
            if len(m.group(0)) >= 5]
    say("Homopolymer runs >= 5: " + (", ".join(f"{b}x{n} @{p}" for b, n, p in runs) or "none"))
    reps = set()
    for k in (12, 10):
        for i in range(L - k + 1):
            w = nt[i:i + k]
            j = nt.find(w, i + 1)
            if j != -1:
                reps.add((k, w, i + 1, j + 1, "direct"))
            r = str(Seq(w).reverse_complement())
            j = nt.find(r)
            if j != -1 and j != i:
                reps.add((k, w, i + 1, j + 1, "inverted"))
        if reps:
            break
    say("Repeats (>=10 nt): " + ("; ".join(f"{t} {w} @{a}&{b}" for k, w, a, b, t in sorted(reps))
                               if reps else "none"))

    text = "\n".join(out)
    (HERE / "report.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
