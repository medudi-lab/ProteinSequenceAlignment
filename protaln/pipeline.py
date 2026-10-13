"""End-to-end pipeline for a TF / partner protein pair across species."""
from .analysis import (compare_to_reference, consensus, disagreement_regions, find_motifs,
                       qc_sequence, species_score)
from .fetch import species_from_header


def best_per_species(records):
    """Keep one sequence per species: the longest."""
    best = {}
    for header, seq in records:
        sp = species_from_header(header)
        if sp not in best or len(seq) > len(best[sp][1]):
            best[sp] = (header, seq)
    return best


def analyse_protein(ref_seq, records, motifs=None):
    """Analyse one protein across species. Returns dict with per-species rows and consensus."""
    ref_seq = ref_seq.rstrip("*")
    rows, projected = [], []
    for sp, (header, seq) in sorted(best_per_species(records).items()):
        seq = seq.rstrip("*")
        qc = qc_sequence(seq)
        cmp_ = compare_to_reference(ref_seq, seq)
        row = {"species": sp, "header": header, "qc": qc, "cmp": cmp_,
               "motifs": find_motifs(seq, motifs), "score": species_score(cmp_, qc)}
        rows.append(row)
        if qc["pass"]:
            projected.append(cmp_["residues"])
    cols = consensus(ref_seq, projected) if projected else []
    return {"rows": rows, "consensus": cols,
            "disagreement_regions": disagreement_regions(cols),
            "reference_motifs": find_motifs(ref_seq, motifs)}


def pair_table(tf_result, partner_result):
    """Join per-species results of the TF and partner; rank by the weaker of the two."""
    tf = {r["species"]: r for r in tf_result["rows"]}
    pr = {r["species"]: r for r in partner_result["rows"]}
    table = []
    for sp in sorted(set(tf) | set(pr)):
        a, b = tf.get(sp), pr.get(sp)
        both = bool(a and b)
        pair_score = min(a["score"], b["score"]) if both else 0.0
        table.append({
            "species": sp,
            "tf_score": a["score"] if a else "", "partner_score": b["score"] if b else "",
            "tf_identity": round(a["cmp"]["identity"], 3) if a else "",
            "partner_identity": round(b["cmp"]["identity"], 3) if b else "",
            "tf_complete": a["cmp"]["completeness"]["complete"] if a else "",
            "partner_complete": b["cmp"]["completeness"]["complete"] if b else "",
            "tf_qc_issues": ",".join(a["qc"]["issues"]) if a else "missing",
            "partner_qc_issues": ",".join(b["qc"]["issues"]) if b else "missing",
            "pair_score": pair_score,
        })
    table.sort(key=lambda r: -r["pair_score"])
    return table
