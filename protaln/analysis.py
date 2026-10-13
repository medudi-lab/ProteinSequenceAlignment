"""QC, motif detection, consensus/disagreement, completeness and scoring."""
import re
from collections import Counter

from .align import align, project, score

VALID = set("ACDEFGHIKLMNPQRSTVWY")
# Common TF-related motifs (regex on the ungapped sequence); users can add their own.
DEFAULT_MOTIFS = {
    "NLS_monopartite": r"K[KR].[KR]",
    "NLS_bipartite": r"[KR]{2}.{10,12}[KR]{3}",
    "C2H2_zinc_finger": r"C.{2,4}C.{12}H.{3,5}H",
    "LxxLL": r"L..LL",
    "bHLH_basic_region": r"[KR].{2}[KR]R.{2}[ER]",
    "N_glycosylation": r"N[^P][ST]",
}


def qc_sequence(seq, min_len=50, max_x_fraction=0.05):
    """Return a dict of QC flags; 'pass' is True when no flag is raised."""
    body = seq.rstrip("*")
    issues = []
    if len(body) < min_len:
        issues.append("too_short")
    if "*" in body:
        issues.append("internal_stop")
    bad = sum(1 for c in body if c not in VALID and c != "*")
    if body and bad / len(body) > max_x_fraction:
        issues.append("ambiguous_residues")
    if not body.startswith("M"):
        issues.append("no_start_Met")
    return {"length": len(body), "ambiguous": bad, "issues": issues,
            "pass": not (set(issues) - {"no_start_Met"})}


def find_motifs(seq, motifs=None):
    """Return {name: [(start_1based, matched_text)]} with overlapping matches."""
    out = {}
    for name, pattern in (motifs or DEFAULT_MOTIFS).items():
        rx = re.compile("(?=(%s))" % pattern)
        out[name] = [(m.start() + 1, m.group(1)) for m in rx.finditer(seq)]
    return out


def completeness(ref_len, residues, insertions):
    """Coverage of the reference and truncation / gap description."""
    covered = [i for i, r in enumerate(residues) if r != "-"]
    if not covered:
        return {"coverage": 0.0, "n_trunc": ref_len, "c_trunc": 0, "internal_gaps": 0,
                "insertions": 0, "complete": False}
    first, last = covered[0], covered[-1]
    internal = sum(1 for r in residues[first:last + 1] if r == "-")
    n_ins = sum(len(v) for k, v in insertions.items() if first <= k < last)
    return {"coverage": len(covered) / ref_len, "n_trunc": first, "c_trunc": ref_len - 1 - last,
            "internal_gaps": internal, "insertions": n_ins,
            "complete": first == 0 and last == ref_len - 1 and internal == 0}


def compare_to_reference(ref, seq):
    """Align seq to ref and compute identity/similarity/completeness/score."""
    ar, aq, raw = align(ref, seq)
    residues, insertions = project(ar, aq)
    pairs = [(a, b) for a, b in zip(residues, ref) if a != "-"]
    aligned = len(pairs)
    ident = sum(1 for a, b in pairs if a == b)
    sim = sum(1 for a, b in pairs if score(a, b) > 0)
    comp = completeness(len(ref), residues, insertions)
    return {"residues": residues, "insertions": insertions, "raw_score": raw,
            "identity": ident / aligned if aligned else 0.0,
            "similarity": sim / aligned if aligned else 0.0,
            "completeness": comp,
            "substitutions": [(i + 1, ref[i], r) for i, r in enumerate(residues)
                              if r not in ("-", ref[i])]}


def consensus(ref, projected, min_fraction=0.5):
    """Column-wise consensus over reference positions.

    projected: list of residue lists (one per species). Returns a list of dicts
    with the consensus residue, its frequency, and whether the column disagrees
    with the reference or is variable.
    """
    cols = []
    for i, refres in enumerate(ref):
        counts = Counter(p[i] for p in projected if p[i] != "-")
        total = sum(counts.values())
        if total:
            top, n = counts.most_common(1)[0]
            freq = n / total
        else:
            top, freq = "-", 0.0
        cols.append({"pos": i + 1, "ref": refres, "consensus": top if freq >= min_fraction else "X",
                     "freq": freq, "n": total,
                     "disagrees_with_ref": total > 0 and top != refres,
                     "variable": freq < min_fraction})
    return cols


def disagreement_regions(cols, min_len=1):
    """Contiguous runs of columns that are variable or disagree with the reference."""
    regions, start = [], None
    flag = [c["variable"] or c["disagrees_with_ref"] for c in cols]
    for i, f in enumerate(flag + [False]):
        if f and start is None:
            start = i
        elif not f and start is not None:
            if i - start >= min_len:
                regions.append((start + 1, i))
            start = None
    return regions


def species_score(cmp_result, qc):
    """Composite 0-100 score: identity, coverage, QC penalty."""
    comp = cmp_result["completeness"]
    s = 100 * (0.5 * cmp_result["identity"] + 0.4 * comp["coverage"] + 0.1 * cmp_result["similarity"])
    s -= 10 * len(qc["issues"])
    if not comp["complete"]:
        s -= 5
    return max(0.0, round(s, 2))
