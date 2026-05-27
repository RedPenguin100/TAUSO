"""Fetch nucleotide records from NCBI (Entrez efetch)."""

import urllib.request

_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"


def fetch_nuccore_fasta(accession: str, timeout: int = 60) -> str:
    """Download a single nuccore record from NCBI efetch and return its FASTA text."""
    with urllib.request.urlopen(_EFETCH_URL.format(acc=accession), timeout=timeout) as resp:
        return resp.read().decode("ascii")
