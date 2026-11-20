"""Offline annotations for DS569k: taxonomic domain (from phylum) and a
protein-family label (normalised protein name). No network access needed."""

import re

import numpy as np
import pandas as pd

ARCHAEA = {"Euryarchaeota", "Thermoproteota", "Candidatus Thermoplasmatota", "Nitrososphaerota",
           "Nanoarchaeota", "Candidatus Lokiarchaeota", "Candidatus Micrarchaeota",
           "Candidatus Heimdallarchaeota", "Candidatus Odinarchaeota", "Candidatus Thorarchaeota"}
EUKARYOTA = {"Chordata", "Streptophyta", "Ascomycota", "Arthropoda", "Nematoda", "Evosea", "Mollusca",
             "Basidiomycota", "Chlorophyta", "Rhodophyta", "Euglenozoa", "Apicomplexa", "Cnidaria",
             "Microsporidia", "Bacillariophyta", "Oomycota", "Ciliophora", "Echinodermata",
             "Platyhelminthes", "Annelida", "Mucoromycota", "Haptophyta", "Cercozoa", "Fornicata",
             "Discosea", "Heterolobosea", "Porifera", "Tardigrada", "Parabasalia", "Chytridiomycota",
             "Nemertea", "Placozoa", "Blastocladiomycota", "Hemichordata", "Brachiopoda",
             "Zoopagomycota", "Tubulinea", "Endomyxa", "Rotifera", "Bryozoa", "Onychophora",
             "Acanthocephala", "Entoprocta", "Chaetognatha", "Priapulida", "Xenacoelomorpha",
             "Ctenophora", "Foraminifera", "Perkinsozoa"}
EUK_CLASSES_NO_PHYLUM = {"Cryptophyceae", "Glaucocystophyceae", "Dinophyceae", "Phaeophyceae",
                         "Raphidophyceae", "Choanoflagellata", "Bigyra", "Olisthodiscophyceae",
                         "Xanthophyceae", "Filasterea", "Chrysophyceae"}
DOMAINS = ["Bacteria", "Eukaryota", "Viruses", "Archaea"]


def domain(phylum, cls, organism):
    if isinstance(phylum, str):
        if phylum.endswith("viricota"):
            return "Viruses"
        if phylum in ARCHAEA:
            return "Archaea"
        if phylum in EUKARYOTA:
            return "Eukaryota"
        return "Bacteria"
    org = str(organism).lower()
    if (isinstance(cls, str) and cls.endswith("viricetes")) or "virus" in org or "phage" in org:
        return "Viruses"
    if cls in EUK_CLASSES_NO_PHYLUM:
        return "Eukaryota"
    return "Unknown"


_UNINFORMATIVE = re.compile(
    r"^(putative |probable )?(uncharacterized|hypothetical|unknown)\b"
    r"|^(putative |probable )?(protein|polypeptide|peptide|orf|gene) [\w.\-/]+$"
    r"|^(putative |probable )?uncharacterised", re.I)
_SUFFIX = re.compile(
    r"(\s*\((fragment|fragments)\))"
    r"|(,\s*(chloroplastic|mitochondrial|peroxisomal|cytoplasmic|plastid|amyloplastic|"
    r"chromoplastic|apicoplast|glyoxysomal|hydrogenosomal|thylakoid|nuclear)\b.*)"
    r"|(\s+(homolog|homologue|isoform|isozyme|paralog)\b.*)", re.I)
_TRAIL = re.compile(r"(\s+([0-9]+[A-Za-z]?|[A-Za-z][0-9]*|[IVX]+|alpha|beta|gamma|delta|"
                    r"[A-Za-z]{1,2}[0-9]+[A-Za-z]?|\d+[-/]\d+|[A-Z]{1,3}\d*-\d+))+$")


def family(name):
    """Normalise a protein name into a family-level label, or None if uninformative."""
    if not isinstance(name, str):
        return None
    n = name.strip()
    n = re.sub(r"^(putative|probable)\s+", "", n, flags=re.I)
    n = _SUFFIX.sub("", n)
    n = re.sub(r"\s*\d+(\.\d+)?\s*kDa\b", "", n)
    n = _TRAIL.sub("", n).strip(" ,-")
    if not n or _UNINFORMATIVE.search(n) or _UNINFORMATIVE.search(name) or len(n) < 4:
        return None
    return n[0].upper() + n[1:]


def annotate(meta):
    out = meta.copy()
    out["domain"] = [domain(p, c, o) for p, c, o in
                     zip(meta.ncbi_taxonomy_phylum, meta.get("ncbi_taxonomy_class",
                         pd.Series([None] * len(meta))), meta.organism_name)]
    out["family_named"] = meta.protein_name.map(family)
    return out

