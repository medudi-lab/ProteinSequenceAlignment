"""Command line interface.

Examples:
  python -m protaln --tf-ref tf_ref.fa --partner-ref partner_ref.fa \\
      --tf-gene SOX2 --partner-gene POU5F1 --outdir results      # mine NCBI (birds)
  python -m protaln --tf-ref tf_ref.fa --partner-ref p_ref.fa \\
      --tf-fasta tf_birds.fa --partner-fasta p_birds.fa --outdir results   # local/VGP FASTA
"""
import argparse
import csv
import os

from .fasta import read_fasta, write_fasta
from .fetch import AVES_TAXID, mine_gene
from .pipeline import analyse_protein, pair_table


def _records(fasta, gene, taxon, api_key, name):
    if fasta:
        return read_fasta(fasta)
    if gene:
        return mine_gene(gene, taxon_id=taxon, api_key=api_key)
    raise SystemExit("Provide --%s-fasta or --%s-gene" % (name, name))


def write_report(outdir, label, result):
    with open(os.path.join(outdir, "%s_per_species.tsv" % label), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "length", "qc_pass", "qc_issues", "identity", "similarity",
                    "coverage", "n_trunc", "c_trunc", "internal_gaps", "complete", "score",
                    "substitutions", "motifs"])
        for r in result["rows"]:
            c, comp = r["cmp"], r["cmp"]["completeness"]
            w.writerow([r["species"], r["qc"]["length"], r["qc"]["pass"], ",".join(r["qc"]["issues"]),
                        "%.3f" % c["identity"], "%.3f" % c["similarity"], "%.3f" % comp["coverage"],
                        comp["n_trunc"], comp["c_trunc"], comp["internal_gaps"], comp["complete"],
                        r["score"], ";".join("%s%d%s" % (a, p, b) for p, a, b in c["substitutions"]),
                        ";".join("%s@%d" % (n, h[0][0]) for n, h in r["motifs"].items() if h)])
    with open(os.path.join(outdir, "%s_consensus.tsv" % label), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "ref", "consensus", "freq", "n", "disagrees_with_ref", "variable"])
        for c in result["consensus"]:
            w.writerow([c["pos"], c["ref"], c["consensus"], "%.3f" % c["freq"], c["n"],
                        c["disagrees_with_ref"], c["variable"]])
    with open(os.path.join(outdir, "%s_disagreement_regions.tsv" % label), "w") as fh:
        fh.write("start\tend\n")
        for s, e in result["disagreement_regions"]:
            fh.write("%d\t%d\n" % (s, e))


def main(argv=None):
    p = argparse.ArgumentParser(prog="protaln", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf-ref", required=True, help="FASTA of the known TF protein")
    p.add_argument("--partner-ref", required=True, help="FASTA of the known partner protein")
    for n in ("tf", "partner"):
        p.add_argument("--%s-gene" % n, help="gene symbol to mine from NCBI")
        p.add_argument("--%s-fasta" % n, help="local FASTA (e.g. VGP assembly translations)")
    p.add_argument("--taxon", type=int, default=AVES_TAXID, help="NCBI taxid (default birds, 8782)")
    p.add_argument("--api-key", help="NCBI API key")
    p.add_argument("--outdir", default="results")
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    results = {}
    for label, ref_path, gene, fasta in (("tf", a.tf_ref, a.tf_gene, a.tf_fasta),
                                         ("partner", a.partner_ref, a.partner_gene, a.partner_fasta)):
        ref = read_fasta(ref_path)[0][1]
        recs = _records(fasta, gene, a.taxon, a.api_key, label)
        with open(os.path.join(a.outdir, "%s_sequences.fa" % label), "w") as fh:
            write_fasta(recs, fh)
        results[label] = analyse_protein(ref, recs)
        write_report(a.outdir, label, results[label])

    table = pair_table(results["tf"], results["partner"])
    with open(os.path.join(a.outdir, "pair_ranking.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0].keys()) if table else ["species"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(table)
    print("Wrote reports for %d species to %s" % (len(table), a.outdir))


if __name__ == "__main__":
    main()
