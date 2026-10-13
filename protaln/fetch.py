"""Retrieve protein sequences from NCBI (Gene/Protein databases, which also host
Vertebrate Genomes Project assemblies' RefSeq/GenBank annotations) via E-utilities.

Vertebrate Genomes Project sequences not yet annotated at NCBI can be supplied
as local FASTA files instead (see the CLI ``--fasta`` option).
"""
import json
import time
import urllib.parse
import urllib.request

from .fasta import parse_fasta

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
AVES_TAXID = 8782


def _get(endpoint, params, api_key=None):
    if api_key:
        params = dict(params, api_key=api_key)
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode()


def build_query(gene, taxon_id=AVES_TAXID, refseq_only=True):
    query = "%s[Gene Name] AND txid%d[Organism:exp]" % (gene, taxon_id)
    if refseq_only:
        query += " AND refseq[filter]"
    return query


def search_proteins(gene, taxon_id=AVES_TAXID, retmax=500, refseq_only=True, api_key=None):
    """Return protein UIDs for a gene symbol within a taxon (default: birds)."""
    data = _get("esearch.fcgi", {"db": "protein", "term": build_query(gene, taxon_id, refseq_only),
                                 "retmax": retmax, "retmode": "json"}, api_key)
    return json.loads(data)["esearchresult"]["idlist"]


def fetch_fasta(ids, batch=100, api_key=None):
    """Return [(header, seq)] for protein UIDs."""
    records = []
    for i in range(0, len(ids), batch):
        text = _get("efetch.fcgi", {"db": "protein", "id": ",".join(ids[i:i + batch]),
                                    "rettype": "fasta", "retmode": "text"}, api_key)
        records.extend(parse_fasta(text))
        time.sleep(0.34)  # NCBI rate limit without an API key
    return records


def mine_gene(gene, **kwargs):
    api_key = kwargs.get("api_key")
    ids = search_proteins(gene, **kwargs)
    return fetch_fasta(ids, api_key=api_key)


def species_from_header(header):
    """Extract '[Species name]' from an NCBI header, else fall back to the first token."""
    if "[" in header and header.rstrip().endswith("]"):
        return header[header.rindex("[") + 1:-1]
    return header.split()[0]
