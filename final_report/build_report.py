"""Final report: is the 288 bp query sequence pathogenic?

Assembles every analysis from the screening exercise into one PDF, with a detailed
methods section per analysis, a results-section caption under every figure, and
explicit limitations.

Usage:  python build_report.py
"""

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
EMB = PROJECT / "esm_embedding" / "results"
COMP = PROJECT / "composition" / "figures"
OUT = HERE / "sequence_pathogenicity_report.pdf"

INK = colors.HexColor("#0b0b0b")
INK_2 = colors.HexColor("#52514e")
INK_MUTED = colors.HexColor("#8a8984")
GRID = colors.HexColor("#e4e3df")
ACCENT = colors.HexColor("#2a78d6")
WARN = colors.HexColor("#e34948")
PAGE_W, PAGE_H = A4
MARGIN = 19 * mm
FRAME_W = PAGE_W - 2 * MARGIN

ss = getSampleStyleSheet()
S = {
    "title": ParagraphStyle("title", parent=ss["Title"], fontName="Helvetica-Bold",
                            fontSize=19, leading=23, textColor=INK, alignment=0,
                            spaceAfter=3),
    "subtitle": ParagraphStyle("subtitle", parent=ss["Normal"], fontName="Helvetica",
                               fontSize=10.5, leading=14, textColor=INK_2, spaceAfter=14),
    "h1": ParagraphStyle("h1", parent=ss["Heading1"], fontName="Helvetica-Bold",
                         fontSize=14, leading=18, textColor=INK, spaceBefore=16,
                         spaceAfter=7),
    "h2": ParagraphStyle("h2", parent=ss["Heading2"], fontName="Helvetica-Bold",
                         fontSize=11, leading=14, textColor=INK, spaceBefore=11,
                         spaceAfter=5),
    "h3": ParagraphStyle("h3", parent=ss["Heading3"], fontName="Helvetica-Oblique",
                         fontSize=9.6, leading=13, textColor=INK_2, spaceBefore=8,
                         spaceAfter=3),
    "body": ParagraphStyle("body", parent=ss["Normal"], fontName="Helvetica",
                           fontSize=9.4, leading=13.4, textColor=INK,
                           alignment=TA_JUSTIFY, spaceAfter=7),
    "caption": ParagraphStyle("caption", parent=ss["Normal"], fontName="Helvetica",
                              fontSize=8.3, leading=11.6, textColor=INK_2,
                              alignment=TA_JUSTIFY, spaceBefore=5, spaceAfter=13),
    "mono": ParagraphStyle("mono", parent=ss["Normal"], fontName="Courier",
                           fontSize=7.4, leading=9.6, textColor=INK_2,
                           spaceBefore=3, spaceAfter=8),
    "ref": ParagraphStyle("ref", parent=ss["Normal"], fontName="Helvetica",
                          fontSize=8.2, leading=11.0, textColor=INK_2,
                          spaceAfter=3.6, leftIndent=14, firstLineIndent=-14),
    "eq": ParagraphStyle("eq", parent=ss["Normal"], fontName="Times-Italic",
                         fontSize=11, leading=17, textColor=INK, alignment=1,
                         spaceBefore=9, spaceAfter=9),
    "eqnote": ParagraphStyle("eqnote", parent=ss["Normal"], fontName="Helvetica",
                             fontSize=8.2, leading=11.4, textColor=INK_2,
                             alignment=TA_JUSTIFY, leftIndent=16, rightIndent=16,
                             spaceAfter=9),
    "key": ParagraphStyle("key", parent=ss["Normal"], fontName="Helvetica",
                          fontSize=9.6, leading=13.6, textColor=INK,
                          alignment=TA_JUSTIFY, spaceAfter=7,
                          leftIndent=9, borderPadding=7, borderWidth=0,
                          backColor=colors.HexColor("#f2f5fa")),
}


REFERENCES = [
    ("altschul1990", "Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990). "
     "Basic local alignment search tool. <i>Journal of Molecular Biology</i> 215(3):403\u2013410."),
    ("camacho2009", "Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, "
     "Madden TL (2009). BLAST+: architecture and applications. "
     "<i>BMC Bioinformatics</i> 10:421."),
    ("karlin1990", "Karlin S, Altschul SF (1990). Methods for assessing the statistical "
     "significance of molecular sequence features by using general scoring schemes. "
     "<i>PNAS</i> 87(6):2264\u20132268."),
    ("henikoff1992", "Henikoff S, Henikoff JG (1992). Amino acid substitution matrices from "
     "protein blocks. <i>PNAS</i> 89(22):10915\u201310919."),
    ("buchfink2021", "Buchfink B, Reuter K, Drost HG (2021). Sensitive protein alignments at "
     "tree-of-life scale using DIAMOND. <i>Nature Methods</i> 18:366\u2013368."),
    ("eddy2011", "Eddy SR (2011). Accelerated profile HMM searches. "
     "<i>PLoS Computational Biology</i> 7(10):e1002195."),
    ("johnson2010", "Johnson LS, Eddy SR, Portugaly E (2010). Hidden Markov model speed "
     "heuristic and iterative HMM search procedure. <i>BMC Bioinformatics</i> 11:431."),
    ("smith1981", "Smith TF, Waterman MS (1981). Identification of common molecular "
     "subsequences. <i>Journal of Molecular Biology</i> 147(1):195\u2013197."),
    ("wood2019", "Wood DE, Lu J, Langmead B (2019). Improved metagenomic analysis with "
     "Kraken 2. <i>Genome Biology</i> 20:257."),
    ("wright1990", "Wright F (1990). The \u2018effective number of codons\u2019 used in a "
     "gene. <i>Gene</i> 87(1):23\u201329."),
    ("sharp1987", "Sharp PM, Li WH (1987). The codon adaptation index \u2014 a measure of "
     "directional synonymous codon usage bias, and its potential applications. "
     "<i>Nucleic Acids Research</i> 15(3):1281\u20131295."),
    ("ikemura1985", "Ikemura T (1985). Codon usage and tRNA content in unicellular and "
     "multicellular organisms. <i>Molecular Biology and Evolution</i> 2(1):13\u201334."),
    ("karlin1997", "Karlin S, Mr\u00e1zek J (1997). Compositional differences within and "
     "between eukaryotic genomes. <i>PNAS</i> 94(19):10227\u201310232."),
    ("nakamura2000", "Nakamura Y, Gojobori T, Ikemura T (2000). Codon usage tabulated from "
     "international DNA sequence databases. <i>Nucleic Acids Research</i> 28(1):292."),
    ("lin2023", "Lin Z, Akin H, Rao R, et al. (2023). Evolutionary-scale prediction of "
     "atomic-level protein structure with a language model. "
     "<i>Science</i> 379(6637):1123\u20131130."),
    ("rives2021", "Rives A, Meier J, Sercu T, et al. (2021). Biological structure and function "
     "emerge from scaling unsupervised learning to 250 million protein sequences. "
     "<i>PNAS</i> 118(15):e2016239118."),
    ("salazar2020", "Salazar J, Liang D, Nguyen TQ, Kirchhoff K (2020). Masked language model "
     "scoring. <i>Proceedings of ACL 2020</i>, 2699\u20132712."),
    ("meier2021", "Meier J, Rao R, Verkuil R, Liu J, Sercu T, Rives A (2021). Language models "
     "enable zero-shot prediction of the effects of mutations on protein function. "
     "<i>NeurIPS</i> 34:29287\u201329303."),
    ("mcinnes2018", "McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold "
     "Approximation and Projection for dimension reduction. <i>arXiv</i>:1802.03426."),
    ("johnson2019", "Johnson J, Douze M, J\u00e9gou H (2019). Billion-scale similarity search "
     "with GPUs. <i>IEEE Transactions on Big Data</i> 7(3):535\u2013547."),
    ("abramson2024", "Abramson J, Adler J, Dunger J, et al. (2024). Accurate structure "
     "prediction of biomolecular interactions with AlphaFold 3. "
     "<i>Nature</i> 630:493\u2013500."),
    ("jumper2021", "Jumper J, Evans R, Pritzel A, et al. (2021). Highly accurate protein "
     "structure prediction with AlphaFold. <i>Nature</i> 596:583\u2013589."),
    ("vankempen2024", "van Kempen M, Kim SS, Tumescheit C, et al. (2024). Fast and accurate "
     "protein structure search with Foldseek. <i>Nature Biotechnology</i> 42:243\u2013246."),
    ("chen2005", "Chen L, Yang J, Yu J, et al. (2005). VFDB: a reference database for bacterial "
     "virulence factors. <i>Nucleic Acids Research</i> 33:D325\u2013D328."),
    ("liu2022", "Liu B, Zheng D, Zhou S, Chen L, Yang J (2022). VFDB 2022: a general "
     "classification scheme for bacterial virulence factors. "
     "<i>Nucleic Acids Research</i> 50(D1):D912\u2013D917."),
    ("uniprot2023", "The UniProt Consortium (2023). UniProt: the Universal Protein "
     "Knowledgebase in 2023. <i>Nucleic Acids Research</i> 51(D1):D523\u2013D531."),
    ("paysan2023", "Paysan-Lafosse T, Blum M, Chuguransky S, et al. (2023). InterPro in 2022. "
     "<i>Nucleic Acids Research</i> 51(D1):D418\u2013D427."),
    ("sayers2022", "Sayers EW, Cavanaugh M, Clark K, et al. (2022). GenBank. "
     "<i>Nucleic Acids Research</i> 50(D1):D161\u2013D164."),
    ("cock2009", "Cock PJA, Antao T, Chang JT, et al. (2009). Biopython: freely available "
     "Python tools for computational molecular biology and bioinformatics. "
     "<i>Bioinformatics</i> 25(11):1422\u20131423."),
    ("balaji2022", "Balaji A, Sapoval M, Seto C, et al. (2022). SeqScreen: accurate and "
     "sensitive functional screening of pathogenic sequences via ensemble learning. "
     "<i>Genome Biology</i> 23:133."),
    ("bird1980", "Bird AP (1980). DNA methylation and the frequency of CpG in animal DNA. "
     "<i>Nucleic Acids Research</i> 8(7):1499\u20131504."),
    ("greenbaum2008", "Greenbaum BD, Levine AJ, Bhanot G, Rabadan R (2008). Patterns of "
     "evolution and host gene mimicry in influenza and other RNA viruses. "
     "<i>PLoS Pathogens</i> 4(6):e1000079."),
]
CITE = {k: i + 1 for i, (k, _) in enumerate(REFERENCES)}


def c(*keys):
    """Superscript citation markers, e.g. c("altschul1990", "camacho2009")."""
    n = ",".join(str(CITE[k]) for k in keys)
    return f"<super><font size='6.6' color='#2a78d6'>[{n}]</font></super>"


def P(t, s="body"):
    return Paragraph(t, S[s])


MAX_FIG_H = PAGE_H - 34 * mm - 34 * mm      # frame height, less room for a caption


def fig(path, caption, width=FRAME_W, keep=True, max_h=MAX_FIG_H):
    """Figure + results-section caption, kept on one page where it fits.

    Scales to the frame in BOTH dimensions: a tall figure (the alignment panels run
    to ~3,500 px) otherwise overflows the page height even at full width.
    """
    p = Path(path)
    if not p.exists():
        return [P(f"[missing figure: {p.name}]", "caption")]
    from PIL import Image as PILImage
    with PILImage.open(p) as im:
        w, h = im.size
    scale = min(width / w, max_h / h)
    iw, ih = w * scale, h * scale
    flow = [Image(str(p), width=iw, height=ih), P(caption, "caption")]
    return [KeepTogether(flow)] if keep and ih < 420 else flow


def table(rows, widths, header=True, size=8.2):
    t = Table(rows, colWidths=widths, hAlign="LEFT")
    st = [
        ("FONT", (0, 0), (-1, -1), "Helvetica", size),
        ("TEXTCOLOR", (0, 0), (-1, -1), INK),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 3.4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.4),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, GRID),
    ]
    if header:
        st += [("FONT", (0, 0), (-1, 0), "Helvetica-Bold", size),
               ("TEXTCOLOR", (0, 0), (-1, 0), INK_2),
               ("LINEBELOW", (0, 0), (-1, 0), 0.9, INK_2)]
    t.setStyle(TableStyle(st))
    return t


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7.4)
    canvas.setFillColor(INK_MUTED)
    canvas.drawString(MARGIN, PAGE_H - 11.5 * mm,
                      "Sequence pathogenicity assessment · 288 bp unknown query")
    # the cover is page 1 of the PDF but carries no number, so the body starts at 1
    canvas.drawRightString(PAGE_W - MARGIN, 11 * mm, f"{doc.page - 1}")
    canvas.setStrokeColor(GRID)
    canvas.setLineWidth(0.4)
    canvas.line(MARGIN, PAGE_H - 13.5 * mm, PAGE_W - MARGIN, PAGE_H - 13.5 * mm)
    canvas.restoreState()


def build():
    doc = BaseDocTemplate(str(OUT), pagesize=A4,
                          leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=18 * mm, bottomMargin=16 * mm,
                          title="Sequence pathogenicity assessment",
                          author="Screening exercise")
    frame = Frame(MARGIN, 16 * mm, FRAME_W, PAGE_H - 34 * mm, id="n")
    # The cover is a full-bleed image with its own margins baked in, so it gets a
    # frame spanning the whole page and no running header/footer.
    cover_frame = Frame(0, 0, PAGE_W, PAGE_H, id="cover",
                        leftPadding=0, rightPadding=0,
                        topPadding=0, bottomPadding=0)
    doc.addPageTemplates([
        PageTemplate(id="cover", frames=[cover_frame]),
        PageTemplate(id="all", frames=[frame], onPage=header_footer),
    ])
    f = []

    # ---------------------------------------------------------------- cover
    cover = HERE / "cover.png"
    if cover.exists():
        f += [Image(str(cover), width=PAGE_W, height=PAGE_H),
              NextPageTemplate("all"), PageBreak()]
    else:
        f += [NextPageTemplate("all")]

    # ---------------------------------------------------------------- title
    f += [P("Is the query sequence pathogenic?", "title"),
          P("A 288 bp unknown nucleotide sequence assessed by homology search, taxonomic "
            "classification, codon-usage analysis, protein language-model embedding, "
            "profile search, sequence-likelihood modelling and structure prediction.", "subtitle")]

    f += [P("<b>Summary of findings</b>", "h2")]
    f += [P(
        "No evidence of pathogenicity was found, and no evidence of natural origin was found "
        "either. The sequence has no detectable homolog in any database searched, is not "
        "classified by either Kraken2 database, shows no similarity to any curated virulence "
        "factor, and its encoded peptides score as statistically indistinguishable from "
        "shuffled sequence under a protein language model. AlphaFold 3 predicts no folded domain "
        "for any ORF, and finds no homolog whatever for the two it models confidently. Its "
        "codon usage is consistent with back-translation from protein using a fixed "
        "one-codon-per-amino-acid table. "
        "<b>The most probable interpretation is a synthetic, designed sequence with no natural "
        "homolog.</b> Because every pathogenicity database consulted is homology-based, the "
        "absence of hits is uninformative about function: these methods can only report "
        "hazard by resemblance, and this sequence resembles nothing in them.", "key")]

    f += [P("The query", "h2"),
          P("<font face='Courier' size='7.4'>"
            "ATGAATCAACATAATACCCAAATCAATAAATTTATCTTTTTTAGTAGCTTCAGAATCATCAGTACCACCCCAATTCAATTT"
            "ATCAACAATAATAGGACCAGCAACAGCAATAGCAACATATTCATCATTAGAAGCAGGAGAATTCAAAACAATACCAACAGC"
            "TTTTTTATTATCAGCAGTATTCCAAGTTTTACAAACAGCACCATTAGAATTACCTTTTTCACAATCACCAGCATGCAAATT"
            "ATCAATACCAGAATTTTTATGAGATCTAATATGAACCAAACCCAT</font>", "body")]

    f += [table([
        ["ORF", "Frame", "Length", "Peptide sequence (N→C)"],
        ["ORF1", "+1", "96 aa", "MNQHNTQINKFIFFSSFRIIST… (partial)"],
        ["ORF2", "+2", "86 aa", "MNIIPKSINLSFLVASESSVPP…"],
        ["ORF3", "+3", "25 aa", "MQENSKQYQQLFYYQQYSKFYKQHH"],
        ["ORF4", "−1", "96 aa", "MGLVHIRSHKNSGIDNLHAGDC… (partial)"],
        ["ORF5", "−2", "47 aa", "WVWFILDLIKILVLIICMLVIV… (partial)"],
    ], [34, 36, 44, FRAME_W - 114])]
    f += [P("Open reading frames identified by NCBI ORFfinder. ORF1 and ORF4 occupy the same "
            "288 bp on opposite strands; both run the full length without an internal stop.",
            "caption")]

    # ------------------------------------------------------- introduction
    # heading bound to the figure: fig() keeps the image with its caption, but a
    # heading left outside that group can still be orphaned at a page foot
    f += [KeepTogether([P("Approach", "h1")] + fig(HERE / "approach.png",
             "<b>Figure 1. The eight analyses and the level of abstraction each one works at.</b> "
             "Methods are grouped by what they actually compare: raw nucleotides, translated "
             "peptides, or learned and structural representations that reach beyond detectable "
             "sequence similarity. The right-hand column states what a positive result from each "
             "method would have established. The grouping matters for interpreting a uniformly "
             "negative outcome: these are not eight repetitions of one test, but eight different "
             "definitions of 'resembles something known', and a sequence can fail the strictest "
             "while still tripping a looser one."))]

    f += [P("Why a sequence is interrogated at several levels at once", "h2")]
    f += [P(
        "Homology, the descent of two sequences from a common ancestor, is not a single "
        "observable but a claim that leaves traces at several levels of biological organisation "
        "simultaneously. Crucially, these traces decay at very different rates. DNA accumulates "
        "silent substitutions rapidly, so two genes sharing an ancestor can become entirely "
        "unrecognisable at the nucleotide level while the proteins they encode remain obviously "
        "similar. Protein sequence in turn diverges faster than the three-dimensional fold it "
        "adopts, which is why structure can preserve a signal of common ancestry long after "
        "sequence alignment has ceased to detect anything. It therefore stands to reason that a "
        "method operating at one level may be entirely blind to a relationship that is plainly "
        "visible at another.", "body")]
    f += [P(
        "This observation motivates the design of the present analysis. Rather than apply a "
        "single method and treat its output as decisive, I interrogated the query at three "
        "levels of abstraction (i.e., the nucleotides themselves, the translated peptides, and "
        "learned representations of protein space), such that each analysis samples a different "
        "and only partly dependent layer of evidence. The value of this design is asymmetric. A "
        "positive result from any single method would be informative on its own; a negative "
        "result from any single method would not, because it would remain consistent with a "
        "relationship visible only at some other level. It is the <i>convergence</i> of "
        "negatives across all three levels that carries weight, and it is for this reason that "
        "the methods below are presented as a ladder rather than a list.", "body")]

    f += [P("Homology search: BLASTn, DIAMOND blastp and jackhmmer", "h2")]
    f += [P(
        "Intuitively, homology search is defined as the problem of determining whether a query "
        "sequence resembles anything in a reference database more closely than chance alone "
        "would permit. The qualifier is essential: any two sequences of sufficient length will "
        "align <i>somewhere</i>, and the task is therefore not to find an alignment but to "
        "establish that a given alignment is improbable under a null model of unrelated "
        "sequence. The three methods used here share this logic, and differ in the alphabet over "
        "which they operate and in the richness of the model they bring to bear.", "body")]
    f += [P(
        "<b>BLASTn</b>" + c("altschul1990", "camacho2009") + " operates on nucleotides. An "
        "exhaustive comparison of the query against every database sequence would be "
        "computationally prohibitive, so BLASTn proceeds heuristically: it first identifies "
        "short exact matches (‘words’, here set to 7 nt to maximise sensitivity on a "
        "short query), and then extends these seeds outwards into longer gapped alignments, "
        "retaining those that score highly. Each surviving alignment is assigned an "
        "<b>E-value</b>, defined as the number of alignments of at least that quality one would "
        "expect to observe by chance in a database of the given size." + c("karlin1990") + " More "
        "concretely, an E-value of 10 indicates that ten such hits are expected from chance "
        "alone, and the result is therefore meaningless; an E-value of 1e−50 indicates that "
        "chance is not a plausible explanation.", "body")]
    f += [P(
        "<b>DIAMOND blastp</b>" + c("buchfink2021") + " applies the same seed-and-extend logic to "
        "translated protein sequence, and is for coding regions considerably the more sensitive "
        "approach. The reason is that comparison at the amino acid level can be chemically "
        "informed. Rather than scoring a match as simply correct or incorrect, DIAMOND employs a "
        "substitution matrix (BLOSUM62)" + c("henikoff1992") + ", in which the replacement of "
        "leucine by isoleucine (chemically similar, and routinely interchanged over evolutionary "
        "time) scores nearly as well as an exact match, whereas the replacement of leucine by "
        "aspartate is penalised. Two proteins may therefore be recognised as related at 25% "
        "identity, at which point the underlying DNA retains no detectable similarity "
        "whatsoever. DIAMOND is an accelerated implementation of this search, which is what "
        "renders a comparison of every ORF against the whole of Swiss-Prot practical.", "body")]
    f += [P(
        "<b>jackhmmer</b>" + c("eddy2011", "johnson2010") + " extends the approach further still, "
        "and is the most sensitive sequence-based method employed here. Pairwise search treats "
        "every position in a protein as equally informative, which is biologically incorrect: "
        "real protein families possess a small number of positions that are absolutely conserved "
        "(often catalytic or structural) alongside many that tolerate almost any residue. "
        "jackhmmer captures this heterogeneity by constructing a <b>profile hidden Markov "
        "model</b>, which records, position by position, which residues are acceptable and how "
        "tolerant each site is to insertion or deletion. It then searches the database with that "
        "profile, admits whatever clears a strict inclusion threshold, rebuilds the profile from "
        "the enlarged set, and iterates.", "body")]
    f += [P(
        "The consequence of this iterative procedure is that homologues may be recruited "
        "<i>through intermediates</i>: a sequence too remote to align directly against the query "
        "may nonetheless align against a profile built from sequences that themselves aligned to "
        "it. jackhmmer therefore detects relationships well beyond the reach of any direct "
        "pairwise comparison, and a search that recruits nothing over five iterations is "
        "correspondingly strong evidence of absence. It establishes not merely that no database "
        "sequence resembles the query, but that none can be reached even through a chain of "
        "intermediate relatives.", "body")]

    f += [P("Taxonomic classification: Kraken2", "h2")]
    f += [P(
        "Where homology search asks what a sequence is related to, taxonomic classification asks "
        "a narrower and more concrete question: from which organism did it come? "
        "<b>Kraken2</b>" + c("wood2019") + " approaches this without recourse to alignment at "
        "all. It decomposes the query into every overlapping subsequence of fixed length (here "
        "35 nt, termed <b>k-mers</b>) and looks each one up in a precomputed hash table mapping "
        "k-mers to the genomes in which they occur. Each k-mer thereby casts a vote for a taxon, "
        "and the query is assigned to the lowest common ancestor consistent with the votes "
        "cast.", "body")]
    f += [P(
        "Because matching is exact, the method is extremely fast but correspondingly brittle: it "
        "cannot recognise a diverged relative, only near-identical sequence. This brittleness is "
        "precisely what makes it useful here, since it constitutes a clean test of a sharply "
        "defined proposition, namely that the query derives from an organism that has already "
        "been sequenced. A result of zero classified k-mers across two databases establishes "
        "that no 35 nt window of the query appears in any catalogued genome.", "body")]

    f += [P("Nucleotide composition and codon usage", "h2")]
    f += [P(
        "Genomes imprint statistical signatures upon their genes, and these signatures persist "
        "independently of what any individual gene does. <b>GC content</b>, the proportion of G "
        "and C bases, varies widely and characteristically between organisms, ranging from "
        "roughly 25% in some <i>Plasmodium</i> species to above 70% in certain actinobacteria. "
        "<b>CpG frequency</b> (a cytosine followed immediately by a guanine) is strongly "
        "suppressed in vertebrate genomes and in many of the viruses that infect "
        "them" + c("bird1980", "greenbaum2008") + ", because methylated cytosine in that context "
        "deaminates to thymine and the dinucleotide is progressively lost; the degree of "
        "suppression is itself informative as to host.", "body")]
    f += [P(
        "<b>Codon usage</b> is the most discriminating of these signals, and it is worth setting "
        "out why. The genetic code is redundant: most amino acids are specified by between two "
        "and six distinct codons, and which of them an organism prefers is shaped both by the "
        "abundance of its tRNA pool and by mutational bias. These preferences are measurable and "
        "species-characteristic" + c("ikemura1985", "karlin1997") + ", which is what allows codon "
        "usage to implicate a source organism even in the complete absence of homology. Two "
        "summary statistics are employed here. The <b>effective number of codons (ENC)</b> "
        "ranges from 61, where all synonymous codons are used with equal frequency, down to 20, "
        "where exactly one codon is used per amino acid" + c("wright1990") + "; real genes "
        "typically fall between 30 and 55. The <b>codon adaptation index (CAI)</b> measures how "
        "closely the codon choices of a gene match the preferred codons of a reference "
        "genome" + c("sharp1987") + ". This analysis proved to be the only one in the present "
        "work to return a positive finding, and, as will become apparent, the one that "
        "ultimately explains the remainder.", "body")]

    f += [P("Protein language-model embedding", "h2")]
    f += [P(
        "Protein language models are trained on tens of millions of natural protein sequences "
        "under a deceptively simple objective: to predict residues that have been masked. "
        "Succeeding at this task requires the model to internalise a great deal of implicit "
        "biology, because which residue belongs at a hidden position depends on its structural "
        "environment, on the residues with which it is in contact, and on what the protein does. "
        "<b>ESM-2</b>" + c("lin2023", "rives2021") + " is such a model, and a sequence passed "
        "through it yields an <b>embedding</b>: a fixed-length vector (1,280 numbers for the "
        "650M-parameter model used here) summarising what the model has inferred about that "
        "protein.", "body")]
    f += [P(
        "The property that makes embeddings useful is that these vectors arrange themselves "
        "meaningfully in the space they occupy. Proteins sharing structure or function come to "
        "lie near one another even where their sequences have diverged past the point at which "
        "alignment detects anything, and nearest-neighbour search in embedding space can "
        "therefore reveal functional relationships invisible to BLAST or even to profile "
        "methods.", "body")]
    f += [P(
        "The weakness of the method matters as much as its strength, and is central to the "
        "interpretation offered in §5. Embeddings are produced by averaging over residue "
        "positions, an operation which necessarily discards order. For short or low-complexity "
        "sequences the resulting average is dominated by amino-acid composition. The "
        "consequence is that a query will always possess nearest neighbours, since the geometry "
        "of the space guarantees as much, whether or not any real relationship exists. Proximity "
        "in embedding space is therefore properly regarded as a hypothesis requiring "
        "independent test, and not as a result in itself.", "body")]

    f += [P("Pairwise alignment of flagged neighbours", "h2")]
    f += [P(
        "Supplying that independent test is the purpose of the pairwise alignment step. Every "
        "embedding neighbour originating from a pathogenic or venomous organism, or itself "
        "matching a curated virulence factor, was realigned against the query using full "
        "<b>Smith-Waterman</b> local alignment" + c("smith1981") + " — that is, the exact "
        "optimal algorithm rather than the heuristic DIAMOND employs, so that weak similarity "
        "cannot be missed through approximation.", "body")]
    f += [P(
        "The crucial element, however, is the manner in which significance is judged. Each "
        "alignment score was compared against a null distribution constructed by shuffling the "
        "query two thousand times and realigning each shuffle in turn. Shuffling preserves "
        "amino-acid composition exactly while destroying order, with the result that a pair "
        "scoring well merely because both sequences are, say, hydrophobic will score equally "
        "well against the shuffles and yield a non-significant p-value. The test thereby "
        "separates genuine homology from compositional coincidence, which is precisely the "
        "confound that mean-pooled embeddings introduce.", "body")]

    f += [P("Pseudo-likelihood: is the sequence protein-like at all?", "h2")]
    f += [P(
        "Every method described thus far asks what the query resembles. The pseudo-likelihood "
        "test asks something different, and arguably more fundamental: whether these peptides "
        "resemble proteins <i>at all</i>, independently of any particular relative. The "
        "procedure masks each residue position in turn, records the probability the model "
        "assigns to the residue actually present, and averages the logarithm of those "
        "probabilities across the sequence." + c("salazar2020", "meier2021") + " The resulting "
        "<b>pseudo-log-likelihood</b> is best understood as a measure of how unsurprising a "
        "sequence is to a model that has learned what natural proteins look like. Real proteins "
        "score highly because their residues sit in contexts evolution has produced repeatedly; "
        "random sequence scores poorly.", "body")]
    f += [P("Formally, for a sequence <i>x</i> = <i>x</i><sub>1</sub>\u2026<i>x</i><sub>L</sub> "
            "of length <i>L</i>, the per-residue pseudo-log-likelihood is", "body")]
    f += [P("PLL(<i>x</i>) &nbsp;=&nbsp; (1/<i>L</i>) &nbsp;\u03a3<sub>i=1..L</sub>&nbsp; "
            "log&nbsp;<i>P</i>(&thinsp;<i>x</i><sub>i</sub> &nbsp;|&nbsp; "
            "<i>x</i><sub>&#8722;i</sub>&thinsp;)", "eq")]
    f += [P(
        "where <i>x</i><sub>&#8722;i</sub> denotes the sequence with position <i>i</i> replaced "
        "by a mask token, so one forward pass is required per position. It is a "
        "<i>pseudo</i>-likelihood because the terms are conditionals computed with the "
        "remainder of the sequence visible, and their sum is therefore not the true joint "
        "log-probability log&nbsp;<i>P</i>(<i>x</i>) \u2014 a quantity a masked language model "
        "never provides directly. Values are negative by construction: a PLL of \u22121.7 "
        "corresponds to the model assigning the observed residue a probability of roughly 0.18 "
        "on average, while \u22122.7 corresponds to about 0.07.", "eqnote")]
    f += [P(
        "Interpreting the resulting score requires care, which is why the test is constructed "
        "around controls. Pseudo-likelihood varies systematically with sequence length, and "
        "every control was therefore length-matched to the ORF against which it is compared. "
        "Four control groups were used: real Swiss-Prot proteins of matched length, serving as "
        "the positive control; those same proteins with their residues shuffled; uniform random "
        "sequence over the twenty amino acids; and — the strictest of the four — the "
        "ORF's own residues shuffled. This last control holds amino-acid composition exactly "
        "constant and varies only the order, and therefore isolates the question of whether the "
        "<i>arrangement</i> of residues is protein-like, rather than merely the ingredients "
        "from which it is composed.", "body")]

    f += [P("Structure prediction", "h2")]
    f += [P(
        "Structure is the most conserved level at which two proteins may resemble one another. "
        "Folds persist across evolutionary distances at which sequence similarity has decayed "
        "entirely, and a predicted structure matching a known toxin or virulence scaffold would "
        "therefore constitute meaningful evidence even in the complete absence of a sequence "
        "hit. <b>AlphaFold 3</b>" + c("abramson2024", "jumper2021") + " predicts "
        "three-dimensional structure from sequence, and reports confidence in two ways.", "body")]
    f += [P(
        "<b>pLDDT</b> (<i>predicted Local Distance Difference Test</i>) scores confidence per "
        "residue on a 0–100 scale, where values above 70 are generally considered reliable "
        "and values below 50 often indicate genuine disorder rather than mere failure. "
        "<b>pTM</b> (<i>predicted Template Modelling score</i>) estimates confidence in the "
        "global fold, running from 0 to 1, where values above roughly 0.5 indicate a "
        "well-defined overall shape. The distinction between them bears directly on the result "
        "reported in §7: pLDDT is a <i>local</i> quantity, scoring a residue's immediate "
        "geometric environment rather than where that piece sits within the whole, whereas pTM "
        "is global. A single well-formed α-helix floating free of any domain may therefore "
        "score above 90 on pLDDT while possessing no tertiary structure whatsoever, and low pTM "
        "is what registers the difference — which is precisely what is observed for ORF3 "
        "and ORF5. AlphaFold additionally reports "
        "<font face='Courier' size='8'>fraction_disordered</font>, the proportion of the chain "
        "it judges disordered, which is how a model may be confident per residue and disordered "
        "overall at one and the same time.", "body")]
    f += [P(
        "One feature of how AlphaFold operates proved as informative as the structures it "
        "returned. The model does not predict from a single sequence in isolation; it first "
        "searches databases for homologues and assembles a multiple-sequence alignment, from "
        "which it extracts constraints encoded in patterns of correlated variation across "
        "relatives. The <b>depth</b> of that alignment, which is to say how many homologues were "
        "found, is reported alongside the result, and therefore constitutes an independent "
        "homology search conducted against databases different from those used in §1.", "body")]

    f += [P("Pathogenicity assessment", "h2")]
    f += [P(
        "Pathogenicity was assessed both directly and indirectly. The direct test searched all "
        "five ORFs against <b>VFDB</b>" + c("chen2005", "liu2022") + ", a curated database of "
        "experimentally characterised bacterial virulence factors — toxins, adhesins, "
        "secretion-system components and the like. The indirect test examined the 250 nearest "
        "embedding neighbours, classifying each according to whether its source organism is a "
        "human pathogen or venomous species and whether the protein itself matches VFDB, and "
        "then subjected every flagged pair to the alignment test described above.", "body")]
    f += [P(
        "Both routes share a limitation to which the conclusion of this report returns at "
        "length. They detect hazard by resemblance to hazards already catalogued, and can "
        "therefore say nothing whatever about a sequence that resembles nothing on record.", "body")]

    # ------------------------------------------------------- 1 homology
    f += [P("1 · Homology search", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "Four independent searches were run. (i) <b>BLASTn</b> against NCBI nt" + c("altschul1990", "sayers2022") + ", with relaxed "
        "settings (word size 7, no low-complexity filtering, expect threshold 1000) to "
        "maximise sensitivity for a short query. (ii) <b>DIAMOND blastp</b> v2.2.9" + c("buchfink2021") + " "
        "(<font face='Courier' size='8'>--very-sensitive</font>) of all five translated ORFs "
        "against all 569,192 Swiss-Prot proteins" + c("uniprot2023") + ", at expect thresholds up to e ≤ 1000. "
        "(iii) <b>DIAMOND blastp</b> of the ORFs against VFDB set A (4,755 curated virulence "
        "factors, downloaded 2026-11-19), e ≤ 10. (iv) <b>jackhmmer</b> (HMMER 3.4)" + c("eddy2011", "johnson2010") + ", five "
        "iterations per ORF against Swiss-Prot, inclusion threshold E &lt; 0.001. Profile "
        "methods were included specifically because iterative HMM search detects remote "
        "homology that pairwise alignment misses.", "body")]
    f += [P("Results", "h3")]
    f += [table([
        ["Search", "Database", "Threshold", "Significant hits"],
        ["BLASTn", "NCBI nt", "E ≤ 1000", "0"],
        ["DIAMOND blastp", "Swiss-Prot (569,192)", "e ≤ 1000", "0  (best e = 9.74, ORF4)"],
        ["DIAMOND blastp", "VFDB set A (4,755)", "e ≤ 10", "0"],
        ["jackhmmer ×5 iter.", "Swiss-Prot (569,192)", "E < 0.001", "0  (best E = 0.7, ORF1)"],
    ], [92, 132, 62, FRAME_W - 286])]
    f += [P(
        "No search returned a significant hit. Of the four, the jackhmmer result is the most "
        "informative, and it is worth stating why. Because no sequence cleared the inclusion "
        "threshold on any iteration, no profile was ever built, and each search therefore "
        "converged immediately upon the query alone. An iterative profile method that cannot "
        "recruit even a single relative from Swiss-Prot constitutes strong evidence that no "
        "detectable homologue exists there — not merely that none aligns directly, but "
        "that none is reachable through intermediates either.", "body")]

    # ------------------------------------------------------- 2 taxonomy
    f += [P("2 · Taxonomic classification", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "Kraken2 v2.1.3 exact 35-mer matching against two databases" + c("wood2019") + ": Standard-8 (archaea, "
        "bacteria, viral, plasmid, human, UniVec) and PlusPF-8 (adding protists and fungi), "
        "both the 2026-10-22 builds. Run in default and <font face='Courier' size='8'>"
        "--report-minimizer-data</font> sensitive modes over the full 288 bp sequence and all "
        "five ORF nucleotide regions.", "body")]
    f += [P("Results", "h3")]
    f += [P(
        "Every query was unclassified in both databases and both modes: <b>0 of 254 k-mers</b> "
        "in the full sequence matched anything, and the same for each ORF region "
        "(100.00% unclassified, 6/6 sequences). A natural sequence of this length from any "
        "sequenced organism would normally share at least some 35-mers with a reference "
        "genome. Sharing none is consistent with a sequence that was never sequenced from an "
        "organism.", "body")]

    # ------------------------------------------------------- 3 composition
    f += [P("3 · Nucleotide composition and codon usage", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "GC content, GC skew, sliding-window GC (60 nt, step 15) and dinucleotide "
        "observed/expected ratios were computed directly. Codon usage was evaluated in both "
        "the +1 and −1 frames against 17 reference organisms from the Kazusa codon-usage" + c("nakamura2000") + " "
        "database, scoring the effective number of codons (ENC, Wright 1990)" + c("wright1990") + ", codon adaptation "
        "index (CAI) and mean synonymous log-likelihood. Two shuffle tests were run: "
        "<b>test A</b> drew 20,000 random sequences of matched base composition to estimate "
        "P(both strands stop-free); <b>test B</b> kept ORF4's protein fixed and redrew "
        "synonymous codons from each reference organism's usage (n = 5,000 each) to estimate "
        "P(opposite strand stop-free) under natural codon usage.", "body")]
    f += [P("Results", "h3")]
    f += [table([
        ["Measure", "Value", "Interpretation"],
        ["GC content", "31.2%", "AT-rich; within natural range"],
        ["CpG obs/exp", "0.00", "zero CpG; vertebrates ≈0.2–0.4, bacteria ≈0.8–1.2"],
        ["ENC, frame −1 (ORF4)", "20.0", "the theoretical minimum — one codon per amino acid"],
        ["GC3, frame −1 (ORF4)", "6.5%", "near-total third-position bias"],
        ["CTA / TCA / TTA usage", "0 of 96", "codons whose reverse complement is a stop"],
        ["Shuffle test A", "P < 0.00005", "both frames stop-free does not happen by chance"],
    ], [112, 68, FRAME_W - 180])]
    f += fig(COMP / "codon_bias.png",
             "<b>Figure 2. Codon usage in the −1 frame is at the theoretical limit of bias.</b> "
             "ORF4 uses exactly one codon for each of the 17 amino acids that have synonymous "
             "alternatives, giving an effective number of codons (ENC) of 20.0 — the floor of "
             "the scale, which no natural coding sequence approaches. Back-translating ORF4's "
             "protein with this 19-entry codon table reproduces the observed 288 bp base for "
             "base. The three codons whose reverse complements are stop codons (CTA, TCA, TTA) "
             "are used zero times in 96 codons, which is what keeps the opposite strand open "
             "and permits ORF1 to run the full length. Shuffle test B shows that preserving "
             "ORF4's protein while drawing codons from natural usage leaves the opposite strand "
             "stop-free only 3.9–10.1% of the time, depending on the reference organism.")
    f += fig(COMP / "composition.png",
             "<b>Figure 3. Base composition is unremarkable; the signal is in codon choice, not "
             "nucleotide content.</b> GC content (31.2%) and sliding-window GC (25–38%) fall "
             "within the range of many AT-rich bacterial and organellar genomes, so composition "
             "alone does not mark the sequence as unusual. The one compositional anomaly is the "
             "complete absence of CpG dinucleotides (observed/expected = 0.00), which follows "
             "mechanically from the fixed codon table rather than from any separate design "
             "choice.")

    # ------------------------------------------------------- 4 embeddings
    f += [P("4 · Protein language-model embedding", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "ORF peptides were embedded with ESM-2 650M" + c("lin2023") + " "
        "(<font face='Courier' size='8'>esm2_t33_650M_UR50D</font>), mean-pooling the final "
        "layer over residue positions with the <font face='Courier' size='8'>&lt;cls&gt;</font> "
        "and <font face='Courier' size='8'>&lt;eos&gt;</font> tokens excluded. This convention "
        "was verified against a precomputed Swiss-Prot embedding set (Zenodo record 22935411) "
        "by re-embedding six proteins "
        "and confirming cosine identity of 1.000000. The supplied set covered 540,881 proteins "
        "but contained no viral entries; the 20,960 missing proteins embedded here (including "
        "12,841 of 17,148 viral entries) were computed with the same convention and merged, "
        "giving 561,841 reference proteins. All vectors were L2-normalised, reduced to 256 "
        "dimensions by PCA, and projected with UMAP (n_neighbors = 15, min_dist = 0.1) on a "
        "FAISS HNSW k-nearest-neighbour graph." + c("mcinnes2018", "johnson2019") + " Neighbour rankings were computed as cosine "
        "similarity in the full 1280-d space, not from UMAP coordinates, which distort local "
        "distance.", "body")]
    f += [P("Results", "h3")]
    f += fig(EMB / "uniprot_650m_umap_4panel.png",
             "<b>Figure 4. The ORFs fall in a region of embedding space populated by small "
             "proteins of every taxonomic domain, not by any single family.</b> "
             "(A) UMAP of all 561,841 Swiss-Prot proteins coloured by taxonomic domain, with "
             "the five ORFs as orange stars; the dashed box marks the region enlarged in C and "
             "D. (B) The same projection coloured by InterPro homologous superfamily" + c("paysan2023") + ", the "
             "SCOP/CATH-derived superfamily classification, for the 20 most abundant "
             "superfamilies (78.6% of proteins carry an assignment). (C) The ORF neighbourhood "
             "enlarged, by domain: bacterial, eukaryotic, archaeal and viral proteins are all "
             "represented, and the ORFs do not sit within any domain-specific cluster. "
             "(D) The same region coloured by superfamily, ranked <i>within</i> this window "
             "rather than globally — the globally dominant superfamilies cover only 3,943 of "
             "the 76,788 proteins here. The superfamilies that do populate this region "
             "(ribosomal proteins, knottin/scorpion-toxin, phospholipase A2, Kunitz domains) "
             "share no function; what they share is small size. The ORFs' nearest neighbours "
             "are therefore best explained by length and amino-acid composition rather than by "
             "homology.")

    f += [P("Nearest-neighbour tables", "h2")]
    f += [P(
        "For each ORF the 50 nearest Swiss-Prot proteins by cosine distance were tabulated with "
        "superfamily, family, abbreviated function, source organism, and two pathogenicity "
        "fields: whether the source organism is a recognised pathogen, and whether the "
        "neighbour protein itself has DIAMOND homology to a VFDB virulence factor "
        "(e ≤ 1e−5, ≥30% identity). A DIAMOND expect value for each ORF↔neighbour pair is "
        "reported alongside the cosine distance.", "body")]
    f += fig(EMB / "neighbour_tables" / "ORF4_top50.png",
             "<b>Figure 5. Embedding proximity without sequence homology: the 50 nearest "
             "neighbours of ORF4.</b> Cosine distances are small (0.0178–0.0311, i.e. cosine "
             "similarity 0.97–0.98), yet the DIAMOND column reads <i>no hit</i> for all 50 "
             "pairs — none align to ORF4 at any expect threshold up to 1000. The neighbour set "
             "is biologically incoherent: archaeal uncharacterised proteins, <i>Bacillus</i> "
             "prophage proteins, phage capsid and tail-fibre proteins, four <i>Chlamydia</i> "
             "virulence-plasmid pGP3 entries, a fungal immunomodulator and an insect "
             "segmentation protein. <b>No coherent protein family is represented</b>, which is "
             "what the compositional explanation predicts and a genuine homology signal would "
             "not. Seventeen neighbours derive from pathogenic organisms "
             "and three carry VFDB homology, but these flags describe the neighbours, not "
             "ORF4, which has no VFDB homology at e ≤ 10. Equivalent tables for ORF1–3 and "
             "ORF5 show the same pattern: 0 of 50 aligned in every case.")

    # ------------------------------------------------------- 5 alignment
    f += [P("5 · Pairwise alignment of flagged neighbours", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "Because DIAMOND is a seed-and-extend heuristic that can miss weak similarity in short "
        "low-complexity sequences, all 92 ORF↔neighbour pairs involving a neighbour from a "
        "pathogenic or venomous source organism, or with VFDB homology, were realigned with "
        "full Smith-Waterman" + c("smith1981") + " (Biopython <font face='Courier' size='8'>PairwiseAligner</font>" + c("cock2009") + ", "
        "local mode, BLOSUM62, gaps −11/−1). Significance was calibrated empirically rather "
        "than analytically: for each pair the ORF was shuffled 2,000 times and realigned, "
        "giving a null distribution of scores from sequences of identical amino-acid "
        "composition. The empirical p-value is the fraction of shuffles scoring at least as "
        "high as the real alignment. This design isolates <i>order</i> from <i>composition</i>: "
        "a pair that aligns well only because both sequences are, for example, hydrophobic will "
        "score equally well against the shuffles and yield a non-significant p.", "body")]
    f += [P("Results", "h3")]
    f += [P(
        "Three of 92 pairs reached nominal p &lt; 0.05 (minimum p = 0.0120). None survives "
        "correction for multiple testing: the Bonferroni threshold for 92 tests is "
        "p &lt; 0.00054, more than twenty-fold below the smallest observed value. The expected "
        "number of false positives at α = 0.05 across 92 tests is 4.6, so observing three is "
        "fewer than chance predicts. Two of the three are the same Vaccinia protein from "
        "different strains, so they are not independent observations. Maximum z-scores per ORF "
        "were +0.79 (ORF1), +1.08 (ORF5), +1.92 (ORF3), +2.66 (ORF4) and +3.13 (ORF2) — several "
        "ORFs align to their flagged neighbours no better than their own shuffled sequence "
        "does.", "body")]
    f += fig(EMB / "alignments" / "ORF2_alignments.png", max_h=PAGE_H - 34 * mm - 62 * mm,
             caption=
             "<b>Figure 6. The alignments contain no conserved motifs.</b> Local Smith-Waterman "
             "alignments of ORF2 against its twelve flagged neighbours, ordered by empirical "
             "p-value; identical residues are tinted and bold, conserved substitutions "
             "(BLOSUM62 &gt; 0) mid-tone, and runs of ≥4 consecutive identities boxed. Genuine "
             "homology produces contiguous conserved blocks corresponding to structural or "
             "catalytic motifs. Across all twelve alignments there is exactly one such run "
             "(<font face='Courier' size='7.5'>PPQF</font>, an N-terminal coincidence). Percent "
             "identity is inversely related to alignment length — the highest values (55%) come "
             "from the shortest alignments (33 residues) while the longest alignment "
             "(73 residues) gives the lowest (37%) — which is the signature of noise rather "
             "than homology. Four neighbours have negative z-scores, meaning ORF2's own "
             "shuffled sequence aligns to them better than ORF2 does, despite cosine distances "
             "placing them among its 50 nearest proteins. Note also that nearly every neighbour "
             "is a small uncharacterised Vaccinia protein of 7.3–9.9 kDa: a size class, not a "
             "protein family.")

    # ------------------------------------------------------- 7 PLL
    f += [P("6 · Are the peptides protein-like? Pseudo-likelihood test", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "Pseudo-log-likelihood (PLL) was computed with ESM-2 650M by masking each residue "
        "position in turn" + c("salazar2020", "lin2023") + " and recording the model's log-probability of the observed residue, "
        "averaged over positions. PLL correlates strongly with sequence length, so every "
        "control was length-matched to the ORF it is compared against. Four control groups of "
        "n = 10 were generated per ORF: (i) <b>real Swiss-Prot proteins</b> of matched length "
        "(±10%), the positive control; (ii) those same proteins with residues shuffled; "
        "(iii) <b>the ORF's own residues shuffled</b> — identical composition, different order, "
        "the strictest control, which asks whether the <i>order</i> is protein-like; and "
        "(iv) uniform random sequence over the 20 amino acids. Percentiles report the fraction "
        "of controls in each group scoring below the ORF.", "body")]
    f += [P("Results", "h3")]
    f += fig(EMB / "pll_controls.png",
             "<b>Figure 7. The ORF peptides are not statistically protein-like.</b> "
             "(A) Each ORF (orange star) against its own length-matched controls. All five ORFs "
             "score between −2.39 and −2.74, far below real Swiss-Prot proteins of the same "
             "length (group mean −1.73) and within the band occupied by the shuffled controls. "
             "(B) Pooled control distributions with the ORFs marked: the ORFs fall inside the "
             "negative-control mass, not the positive-control mass. (C) PLL against length, "
             "confirming that the separation is not a length artefact — the ORFs sit below the "
             "real-protein cloud at every length. Against their own composition-matched "
             "shuffles, ORF1 and ORF5 sit at the 50th percentile, meaning <b>their own residues "
             "in random order are exactly as protein-like as the real sequence</b>; ORF3 and "
             "ORF4 reach the 100th percentile but by a visibly small margin. Against real "
             "proteins, every ORF falls in the bottom 40% and ORF2 and ORF4 in the bottom 10%. "
             "A model trained on 65 million natural sequences does not recognise these peptides "
             "as the kind of sequence evolution produces.")

    # ------------------------------------------------------- 8 structure
    f += [P("7 · Structure prediction (AlphaFold 3)", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "All five ORFs were submitted to the AlphaFold 3 server" + c("abramson2024") + ", which returns five seeds per "
        "job together with the multiple-sequence alignments (MSAs) and structural templates it "
        "assembled. Per-residue confidence (pLDDT) was read from the B-factor column of each "
        "mmCIF model and averaged per residue over atoms; pTM, ipTM and "
        "<font face='Courier' size='8'>fraction_disordered</font> were taken from the server's "
        "summary files. MSA depth was counted directly from the returned "
        "<font face='Courier' size='8'>.a3m</font> files — this is the number of homologous "
        "sequences AlphaFold itself was able to find, and so is an independent homology search "
        "run against AlphaFold's own databases rather than the ones used in §1. Secondary "
        "structure was quantified geometrically as the fraction of i→i+4 Cα distances falling "
        "in the α-helical range (5.0–7.5 Å), and chain compactness as the end-to-end distance "
        "per residue, where ≈1.5 Å/residue indicates a fully extended rod rather than a folded "
        "chain. Foldseek structural search was not run, for the reason given in the Results.",
        "body")]
    f += [P(
        "One submission artefact should be recorded. The ORF2 job was submitted as a two-chain "
        "complex in which chain A repeats ORF1's sequence and chain B is the true ORF2; it is "
        "consequently the only job carrying an ipTM (0.08). All ORF2 values reported here are "
        "taken from chain B. The duplicated chain A reproduces ORF1's own result closely (mean "
        "pLDDT 32.3 in the complex against 41.6 alone), so nothing is lost, but ORF2 has not "
        "been modelled as an isolated monomer.", "body")]
    f += [P("Results", "h3")]
    f += fig(EMB / "af3_structures.png",
             "<b>Figure 8. AlphaFold finds no homologs for the two ORFs it is confident about, "
             "and both are single helices rather than folds.</b> "
             "(A) Per-residue pLDDT with all five seeds overlaid; the seeds are nearly "
             "superimposed (mean-pLDDT spread ≤ 2.7 for every ORF), so these are stable "
             "predictions rather than noise. ORF1, ORF2 and ORF4 sit in the 'very low' band "
             "along their entire length, with not a single residue above pLDDT 70 in any of the "
             "three. ORF3 and ORF5 are confident throughout — mean pLDDT 83.6 and 89.5, with "
             "24/25 and 47/47 residues above 70. "
             "(B) The crux of the panel: those two confident models are precisely the two for "
             "which AlphaFold assembled an MSA of depth 1, meaning it found <i>no</i> homologous "
             "sequence anywhere in its databases and the alignment contains only the query "
             "itself. The confidence is therefore single-sequence prediction, not evolutionary "
             "inference, and it independently reproduces the negative homology results of §1 "
             "against a different set of databases. The three ORFs that did recruit some "
             "alignment depth (6–32 sequences) are the three AlphaFold cannot fold. "
             "(C) What the confident models actually are: ORF3 is 90% helical and ORF5 100% "
             "helical, both extended at ≈1.5 Å per residue — straight rods with no tertiary "
             "contacts. AlphaFold's own "
             "<font face='Courier' size='8'>fraction_disordered</font> is 1.00 for both, which "
             "is not a contradiction of their high pLDDT but its explanation: the model is "
             "confident about local backbone geometry and there is no fold to be confident "
             "about. ORF5's sequence (WVWFILDLIKILVLIICMLVIVKKVILMVLFVKLGILLIIKKLLVLF) is a "
             "canonical hydrophobic transmembrane helix. "
             "(D) All metrics tabulated. No ORF produces a folded domain: three are predicted "
             "disordered at low confidence, and two are predicted, confidently, to be lone "
             "helices.")
    f += [P(
        "Foldseek was not run, and on this evidence it would not be informative." + c("vankempen2024") + " Foldseek "
        "matches tertiary folds, and the only two models with usable confidence have no "
        "tertiary structure to match — a single helix is the most common local motif in the PDB "
        "and will return large numbers of high-scoring, functionally unrelated hits from every "
        "structural class. The three models that do have predicted tertiary arrangement are at "
        "pLDDT ≈ 30–42, well below the threshold at which predicted coordinates are reliable "
        "enough to search with, so any Foldseek hit would be a match to noise. This is a real "
        "gap in coverage rather than a closed question, and it is carried into §9.", "body")]
    f += [P(
        "The structural result converges with every other line of evidence and adds one new "
        "observation. The convergence: AlphaFold's inability to find a single homolog for ORF3 "
        "and ORF5 is a fourth independent negative homology search, run against different "
        "databases from BLASTn, DIAMOND and jackhmmer. The new observation: the two peptides "
        "that are confidently modelled are a hydrophobic transmembrane-like helix (ORF5) and a "
        "short polar helix (ORF3), which is a structural character consistent with the "
        "compositional extremes already noted in §3 and §4, and is the kind of sequence that "
        "can be written down without resembling anything natural.", "body")]

    # ------------------------------------------------------- 9 pathogenicity
    f += [P("8 · Pathogenicity assessment", "h1")]
    f += [P("Methods", "h3")]
    f += [P(
        "Pathogenicity was approached from two directions. <b>Direct:</b> all five ORFs were "
        "searched against VFDB set A (4,755 curated virulence factors) with DIAMOND " + c("liu2022") + ""
        "<font face='Courier' size='8'>--very-sensitive</font> at e ≤ 10. <b>Indirect:</b> for "
        "each of the 250 nearest embedding neighbours, the source organism was classified as "
        "pathogenic, venomous or neither, and the neighbour protein itself was searched against "
        "VFDB (e ≤ 1e−5, ≥30% identity); all flagged pairs were then realigned as described in "
        "§5. SeqScreen" + c("balaji2022") + " was evaluated for FunSOC assignment but not run — no hosted service "
        "exists, and local installation requires ~25 dependencies, several without Apple "
        "Silicon builds, against 26 GB of free disk.", "body")]
    f += [P("Results", "h3")]
    f += [table([
        ["Test", "Result"],
        ["ORFs vs VFDB (DIAMOND, e ≤ 10)", "0 hits for all five ORFs"],
        ["Neighbours from pathogenic organisms", "85 of 250 (ORF5 highest at 28/50)"],
        ["Neighbours with VFDB homology", "4 of 249 distinct proteins"],
        ["ORF↔flagged-neighbour alignments", "0 of 92 significant after correction"],
    ], [196, FRAME_W - 196])]
    f += [Spacer(1, 9)]
    f += [P(
        "<b>No ORF has detectable similarity to any known virulence factor.</b> Pathogen-source "
        "and VFDB flags appear among the embedding neighbours, but no ORF aligns to any flagged "
        "neighbour above chance, so no functional inference transfers. The correct reading of "
        "the neighbour tables is that embedding proximity to a virulence protein is not "
        "evidence of a relationship to it.", "body")]
    f += [P(
        "The structural evidence of §7 does not alter this assessment, and sharpens one part of "
        "it. Structure is conserved beyond the point where sequence similarity becomes "
        "undetectable, so a predicted fold matching a toxin or virulence scaffold would have "
        "been the strongest available evidence of hazard; no such fold exists to test, because "
        "no ORF is predicted to adopt a folded domain at usable confidence. The one structural "
        "feature that is predicted confidently — ORF5's hydrophobic transmembrane-like helix — "
        "is a membrane-association motif rather than a virulence determinant. It is shared by "
        "pore-forming toxins and by the great majority of ordinary membrane proteins alike, and "
        "carries no functional specificity on its own.", "body")]

    # ------------------------------------------------------- limitations
    f += [PageBreak(), P("9 · Limitations", "h1")]
    for h, t in [
        ("Short sequences limit every sequence-based method",
         "ORF3 is 25 aa and ORF5 is 47 aa. Below roughly 50–100 residues, local alignment "
         "cannot reach statistical significance even when a genuine relative exists, because "
         "there is insufficient signal to clear the noise threshold in a database of this size. "
         "For these two ORFs, 'no detectable homolog' is a statement about detection limits, "
         "not a demonstration of absence. The conclusion is materially stronger for ORF1, ORF2 "
         "and ORF4 at 86–96 aa."),
        ("Absence of evidence is not evidence of safety",
         "Every pathogenicity resource consulted — VFDB, and SeqScreen's FunSOC framework had "
         "it been run — assigns hazard by homology to curated examples. A sequence with no "
         "homolog is unclassifiable by construction, and would return clean whether it were "
         "harmless or a de novo designed toxin. These results establish that the sequence does "
         "not resemble known hazards; they cannot establish that it is not hazardous."),
        ("The pseudo-likelihood test used n = 10 per group",
         "Percentile resolution is therefore ±10%, sufficient to show clear group separation "
         "but not to place an individual ORF precisely. The most consequential single number — "
         "ORF1 and ORF5 sitting at the 50th percentile against their own shuffles — would "
         "warrant confirmation at higher n before being relied on quantitatively. The original "
         "design specified n = 100; this was reduced for runtime."),
        ("The reference embedding set was incomplete and partly recomputed",
         "The precomputed Swiss-Prot embeddings excluded all 17,148 viral proteins and ~11,000 "
         "others. Of the 28,311 missing, 20,960 were embedded here and merged; the run was "
         "stopped early, leaving 4,307 viral and ~3,000 other proteins (all ≥515 aa) absent. "
         "Since the ORFs are 25–96 aa, the omitted long proteins are unlikely to affect "
         "neighbour rankings, but the reference set is not complete Swiss-Prot."),
        ("Pathogen-organism classification was hand-curated and coarse",
         "The pathogen flag was assigned over the 54 genera and 54 viral species that actually "
         "appear among the neighbours, not from a comprehensive authority. It operates at genus "
         "level, so <i>Bacillus</i> covers both <i>B. subtilis</i> and <i>B. anthracis</i>. It "
         "is adequate for flagging rows in these tables and should not be reused as a general "
         "classifier."),
        ("Mean-pooled embeddings discard sequence order",
         "The embedding analysis in §4 averages over residue positions, which for short "
         "low-complexity sequences is dominated by amino-acid composition. This is itself a "
         "finding — it explains the incoherent neighbour sets — but it means the embedding "
         "results should not be read as evidence about homology in either direction. "
         "Order-aware embedding methods (pLM-BLAST, PROST) were not applied."),
        ("Structures were predicted but never searched",
         "AlphaFold 3 models were obtained (§7), but no structural similarity search was run "
         "against them. Foldseek was judged uninformative on these particular models — the two "
         "confident ones are single helices with no tertiary structure to match, and the three "
         "with tertiary arrangement are at pLDDT ≈ 30–42, too low to search with. That "
         "reasoning is sound for whole-chain fold matching, but it is not a proof that nothing "
         "would be found: Foldseek in TM-align mode, or a search restricted to the ORF5 helix "
         "against transmembrane toxin structures, might still return something interpretable. "
         "Structure remains the most sensitive signal for remote relationships, and this is the "
         "largest untested avenue in the analysis."),
        ("ORF2 was modelled only inside a two-chain job",
         "The ORF2 AlphaFold submission paired it with a duplicate of ORF1's sequence, so ORF2 "
         "has no isolated-monomer prediction. Its chain-B confidence is very low throughout "
         "(mean pLDDT 28.7, 0/86 residues above 70), and a monomer run would be very unlikely "
         "to change that, but the comparison with the other ORFs is not strictly like-for-like. "
         "The ipTM of 0.08 for that job describes an interface between two peptides that have "
         "no reason to interact and should not be read as a result."),
    ]:
        f += [P(f"<b>{h}</b>", "h2"), P(t, "body")]

    # ------------------------------------------------------- conclusion
    f += [P("10 · Conclusion", "h1")]
    f += [P(
        "<b>The question posed was whether this sequence is pathogenic. The answer supported by "
        "the evidence is that it shows no signature of pathogenicity, and no signature of "
        "natural origin either.</b>", "key")]
    f += [P(
        "Eight independent lines of evidence converge upon this conclusion. The sequence "
        "returns no BLASTn hit, no Kraken2 classification in either of two databases (0 of 254 "
        "k-mers), no DIAMOND alignment to Swiss-Prot or to VFDB, and no jackhmmer profile "
        "recruitment over five iterations. Its nearest neighbours in protein language-model "
        "embedding space are explicable by length and composition alone, and none of them "
        "aligns to it. Its encoded peptides score as no more protein-like than their own "
        "shuffled residues. AlphaFold 3 independently finds no homologue for the two peptides "
        "it is able to model confidently (an MSA of depth 1, which is to say the query alone), "
        "and predicts no folded domain for any of the five.", "body")]
    f += [P(
        "Against this uniformly negative background, the codon analysis supplies the one "
        "positive finding, and it proves to be the decisive one. ORF4 employs exactly one codon "
        "for every multi-codon amino acid (ENC = 20.0, the theoretical floor); back-translating "
        "its protein with that 19-entry table reproduces the 288 bp exactly; and the three "
        "codons whose reverse complements are stop codons are avoided entirely, which is "
        "precisely the constraint that permits both strands to remain open across the full "
        "length. No natural coding sequence behaves in this manner. <b>The sequence was "
        "designed, by back-translation from a protein sequence, under a constraint that kept "
        "both reading frames open.</b>", "body")]
    f += [P(
        "It is worth emphasising that this finding explains the pattern of results rather than "
        "merely restating it. Nucleotide-level methods fail because the codon table was chosen "
        "rather than inherited, with the consequence that no natural k-mer history exists "
        "against which to match. Protein-level methods fail because the encoded peptides have "
        "no natural relatives to recover. Embedding methods return neighbours because the "
        "geometry of the space obliges them to return something, and in the absence of genuine "
        "homology what they return is sequences of similar size and composition.", "body")]
    f += [P(
        "For biosecurity purposes, the operative conclusion is therefore a limitation rather "
        "than a clearance, and the distinction is not a semantic one. Every screening resource "
        "applied here detects hazard by resemblance to catalogued examples, and this sequence "
        "resembles nothing catalogued. It follows that a designed sequence possessing no "
        "natural homologue is invisible to homology-based screening whether it is benign or "
        "otherwise — the screen cannot distinguish the two cases, because both present as "
        "silence. Structure prediction, added here precisely in the hope of circumventing this "
        "problem, narrows it less than one might wish: AlphaFold returns no fold to compare "
        "against catalogued hazards, and the single confidently predicted feature (a "
        "transmembrane-like helix in ORF5) is shared by pore-forming toxins and entirely "
        "ordinary membrane proteins alike. Establishing what, if anything, these peptides "
        "actually do would require evidence this analysis cannot supply: a structural "
        "similarity search that the present models cannot support, and ultimately experimental "
        "characterisation.", "body")]

    # ------------------------------------------------------- appendix
    f += [P("Appendix · Software and data", "h1")]
    f += [table([
        ["Component", "Version / source"],
        ["ESM-2", "esm2_t33_650M_UR50D (650M), esm2_t6_8M_UR50D (8M)"],
        ["DIAMOND", "2.2.9"],
        ["HMMER (jackhmmer)", "3.4"],
        ["Kraken2 databases", "Standard-8 and PlusPF-8, 2026-10-22 builds"],
        ["VFDB", "set A core dataset, retrieved 2026-11-19 (4,755 proteins)"],
        ["InterPro", "entry.list, 3,519 homologous-superfamily entries"],
        ["Swiss-Prot reference", "DS569k (569,192 proteins); 561,841 embedded \u2014 "
                                 "huggingface.co/datasets/donnyb/DS569k"],
        ["Precomputed embeddings", "zenodo.org/records/22935411"],
        ["UMAP / FAISS", "umap-learn with FAISS HNSW precomputed kNN"],
        ["Alignment", "Biopython PairwiseAligner, BLOSUM62, gaps −11/−1"],
        ["AlphaFold", "AlphaFold 3 server, 5 seeds per ORF"],
    ], [128, FRAME_W - 128])]
    f += [P(
        "Analysis scripts: <font face='Courier' size='8'>composition/composition_codon.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/uniprot_650m_map.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/neighbour_tables.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/pairwise_flagged.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/alignment_figures.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/pseudolikelihood.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/pll_figure.py</font>, "
        "<font face='Courier' size='8'>esm_embedding/af3_figure.py</font>. "
        "Full per-ORF neighbour tables and alignment figures accompany this report in "
        "<font face='Courier' size='8'>esm_embedding/results/</font>.", "caption")]

    # ------------------------------------------------------- references
    f += [P("References", "h1")]
    f += [P("Numbered as cited in the text.", "caption")]
    for i, (_, text) in enumerate(REFERENCES, 1):
        f += [Paragraph(f"{i}.&nbsp;&nbsp;{text}", S["ref"])]

    # ------------------------------------------------------- AI usage
    f += [KeepTogether([
        P("Statement on AI usage", "h1"),
        P("Claude Code was used throughout this project: to write and debug the analysis "
          "scripts, run the searches and models, generate the figures, and draft this "
          "report. The direction of the project, the experimental design and the key "
          "judgements were my own.", "body"),
    ])]

    doc.build(f)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    build()
