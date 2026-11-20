# Is this sequence pathogenic?

A biosecurity screening exercise: an unknown **288 bp nucleotide sequence** with no
supplied provenance, assessed for evidence of pathogenicity using eight independent
computational approaches.

**[→ Read the full report (PDF, 18 pages)](final_report/sequence_pathogenicity_report.pdf)**

---

## The finding

**No evidence of pathogenicity. No evidence of natural origin either.**

The sequence has no detectable homolog in any database searched, is unclassified by
Kraken2, shows no similarity to any curated virulence factor, and its encoded
peptides score as statistically indistinguishable from shuffled sequence under a
protein language model. AlphaFold 3 predicts no folded domain for any ORF, and finds
no homolog at all for the two it models confidently.

The one positive finding is in codon usage. The −1 frame has an effective number of
codons (ENC) of **20.0** — the theoretical floor, meaning exactly one codon is used
per amino acid. Back-translating ORF4's peptide with that 19-entry table reproduces
the 288 bp exactly, and the three codons whose reverse complements are stop codons
(CTA, TCA, TTA) are absent entirely, which is what keeps both strands open across the
full length. No natural coding sequence behaves this way.

**The most probable interpretation is a designed sequence, back-translated from
protein, with no natural homolog.**

The operative biosecurity point is a limitation rather than a clearance: every
resource used here detects hazard by *resemblance* to catalogued examples, so a
designed sequence with no natural homolog is invisible to all of them — whether it is
benign or not. The absence of hits is uninformative about function.

## Results at a glance

| Analysis | Tool | Result |
|---|---|---|
| Nucleotide homology | BLASTn (NCBI nt, E ≤ 1000) | 0 hits |
| Protein homology | DIAMOND blastp (Swiss-Prot, 569,192) | 0 hits (best e = 9.74) |
| Profile search | jackhmmer ×5 iterations | 0 recruited; no profile ever built |
| Taxonomic classification | Kraken2 (Standard-8, PlusPF-8) | 0 / 254 k-mers classified |
| Virulence factors | DIAMOND vs VFDB set A | 0 hits |
| Embedding neighbours | ESM-2 650M, 561,841 proteins | neighbours explained by composition |
| Flagged-pair alignment | Smith-Waterman, 2,000 shuffles | 0 / 92 significant after Bonferroni |
| Protein-likeness | ESM-2 pseudo-likelihood | ORFs within negative-control band |
| Structure | AlphaFold 3, 5 seeds per ORF | no folded domain; MSA depth 1 for ORF3/ORF5 |
| Codon usage | ENC, CAI, back-translation | **ENC = 20.0; exact back-translation** |

## Repository layout

```
screening-exercise/
├── final_report/                 # the report and the code that builds it
│   ├── sequence_pathogenicity_report.pdf
│   ├── build_report.py           # assembles the PDF; holds all prose + bibliography
│   ├── cover_figure.py           # infographic cover page
│   └── approach_figure.py        # Figure 1, the methods overview
│
├── composition/                  # nucleotide composition and codon usage
│   ├── composition_codon.py      # GC, CpG, ENC, CAI, back-translation test
│   └── figures/                  # Figures 2–3
│
├── esm_embedding/                # protein language-model work
│   ├── uniprot_650m_map.py       # embeds Swiss-Prot, builds the UMAP (Figure 4)
│   ├── embed_missing.py          # embeds proteins absent from the supplied set
│   ├── neighbour_tables.py       # top-50 neighbour tables (Figure 5)
│   ├── pairwise_flagged.py       # Smith-Waterman + shuffle-calibrated p-values
│   ├── alignment_figures.py      # renders the alignments (Figure 6)
│   ├── pseudolikelihood.py       # masked-position PLL with 4 control groups
│   ├── pll_figure.py             # Figure 7
│   ├── af3_figure.py             # Figure 8, from the AlphaFold outputs
│   └── results/                  # result tables and figures
│
├── kraken/                       # Kraken2 reports (both databases, both modes)
├── alphafold_orfs/               # AlphaFold 3 server outputs, 5 seeds per ORF
├── blastn-result-online.xml      # BLASTn result
├── screening-exercise-orfs.fasta # the five ORFs
└── sequence.fasta                # the 288 bp query
```

## The query

```
>sequence
ATGAATCAACATAATACCCAAATCAATAAATTTATCTTTTTTAGTAGCTTCAGAATCATCAGTACCACCCCAATTCAATTT
ATCAACAATAATAGGACCAGCAACAGCAATAGCAACATATTCATCATTAGAAGCAGGAGAATTCAAAACAATACCAACAGC
TTTTTTATTATCAGCAGTATTCCAAGTTTTACAAACAGCACCATTAGAATTACCTTTTTCACAATCACCAGCATGCAAATT
ATCAATACCAGAATTTTTATGAGATCTAATATGAACCAAACCCAT
```

Five ORFs were identified by NCBI ORFfinder. ORF1 (+1, 96 aa) and ORF4 (−1, 96 aa)
occupy the same 288 bp on opposite strands, both running the full length without an
internal stop — which is the feature the codon analysis ultimately explains.

## Reproducing the analysis

Large inputs are **not** committed (see `.gitignore`); the repository keeps code,
small result tables and figures, and weighs ~55 MB rather than ~30 GB. Everything
excluded is either a third-party download or an intermediate the scripts rebuild.

| Excluded | Size | How to obtain |
|---|---|---|
| `DS569k.parquet` | 419 MB | [huggingface.co/datasets/donnyb/DS569k](https://huggingface.co/datasets/donnyb/DS569k) |
| `kraken/kraken_db/` | 22 GB | Standard-8 and PlusPF-8 from the Kraken2 index page |
| `esm_embedding/model_650m/` | 2.4 GB | `facebook/esm2_t33_650M_UR50D`, fetched from HuggingFace |
| `esm_embedding/downloaded_embeddings/` | 2.8 GB | [zenodo.org/records/22935411](https://zenodo.org/records/22935411) |
| `esm_embedding/diamond/` | 449 MB | built by `diamond makedb` from Swiss-Prot and VFDB |

```bash
python -m venv .venv && source .venv/bin/activate
pip install torch transformers numpy pandas pyarrow biopython \
            matplotlib umap-learn faiss-cpu reportlab pillow
```

Then, in order:

```bash
python composition/composition_codon.py          # composition, codon usage
python esm_embedding/uniprot_650m_map.py all     # embed, UMAP, plot
python esm_embedding/neighbour_tables.py         # top-50 neighbour tables
python esm_embedding/pairwise_flagged.py         # alignments + empirical p-values
python esm_embedding/alignment_figures.py        # render the alignments
python esm_embedding/pseudolikelihood.py --n 10  # PLL against matched controls
python esm_embedding/pll_figure.py
python esm_embedding/af3_figure.py               # needs alphafold_orfs/extracted/
python final_report/cover_figure.py
python final_report/approach_figure.py
python final_report/build_report.py              # assembles the PDF
```

`af3_figure.py` reads the AlphaFold outputs unpacked into
`alphafold_orfs/extracted/<orf>/`; unzip the five archives there first.

### Platform note

On macOS, `torch` and `faiss` cannot be imported in the same process — the OpenMP
runtimes clash and segfault. The pipeline is split into separate stages for this
reason, and they must be run as separate processes rather than combined.

## Known limitations

The report has a full limitations section; the three that most affect
interpretation:

- **Short sequences.** ORF3 is 25 aa and ORF5 is 47 aa. Below roughly 50–100
  residues, local alignment cannot reach statistical significance even when a
  genuine relative exists.
- **Structures were predicted but never searched.** Foldseek was judged
  uninformative on these particular models — the two confident ones are single
  helices with no tertiary structure to match, and the other three are at
  pLDDT ≈ 30–42 — but that is a judgement, not a proof that nothing would be found.
  This is the largest untested avenue.
- **The PLL test ran at n = 10** per control group rather than the n = 100
  originally specified. The ORF1/ORF5 50th-percentile result should be confirmed at
  higher n before being relied on quantitatively.

## Statement on AI usage

Claude Code was used throughout this project: to write and debug the analysis
scripts, run the searches and models, generate the figures, and draft the report. The
direction of the project, the experimental design and the key judgements were my own.
