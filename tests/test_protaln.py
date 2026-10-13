import os
import tempfile
import unittest

from protaln.align import align, project
from protaln.analysis import (compare_to_reference, consensus, disagreement_regions,
                              find_motifs, qc_sequence)
from protaln.cli import main
from protaln.fasta import parse_fasta
from protaln.pipeline import analyse_protein

REF = "MSKRKAPQETLNGGITDMLTELANFEKNVSQAIHKYNAGDMSSSNF"


class Tests(unittest.TestCase):
    def test_fasta(self):
        self.assertEqual(parse_fasta(">a x\nAC\nGT\n>b\nM"), [("a x", "ACGT"), ("b", "M")])

    def test_truncated_and_substitution(self):
        q = REF[5:30].replace("E", "D", 1)
        c = compare_to_reference(REF, q)
        self.assertEqual(c["completeness"]["n_trunc"], 5)
        self.assertEqual(c["completeness"]["c_trunc"], len(REF) - 30)
        self.assertFalse(c["completeness"]["complete"])
        self.assertEqual(len(c["substitutions"]), 1)

    def test_internal_deletion(self):
        c = compare_to_reference(REF, REF[:20] + REF[24:])
        self.assertEqual(c["completeness"]["internal_gaps"], 4)

    def test_insertion_projection(self):
        ar, aq, _ = align(REF, REF[:20] + "WWW" + REF[20:])
        res, ins = project(ar, aq)
        self.assertEqual("".join(res), REF)
        self.assertEqual(sum(len(v) for v in ins.values()), 3)

    def test_qc(self):
        self.assertTrue(qc_sequence(REF)["pass"] or "too_short" in qc_sequence(REF)["issues"])
        self.assertIn("internal_stop", qc_sequence("M" + "A" * 60 + "*AAA")["issues"])
        self.assertIn("ambiguous_residues", qc_sequence("M" + "X" * 60)["issues"])

    def test_motifs(self):
        self.assertEqual(find_motifs("AAKRKRAA", {"nls": "K[KR].[KR]"})["nls"][0][0], 3)

    def test_consensus(self):
        proj = [list("ABC"), list("ABC"), list("AXC")]
        cols = consensus("ABD", proj)
        self.assertTrue(cols[2]["disagrees_with_ref"])
        self.assertEqual(disagreement_regions(cols), [(3, 3)])

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            def w(name, text):
                path = os.path.join(d, name)
                open(path, "w").write(text)
                return path
            tf = w("tf.fa", ">ref\n%s\n" % REF)
            birds = w("b.fa", ">1 x [Gallus gallus]\n%s\n>2 x [Taeniopygia guttata]\n%s\n"
                      % (REF, REF[:30]))
            main(["--tf-ref", tf, "--partner-ref", tf, "--tf-fasta", birds,
                  "--partner-fasta", birds, "--outdir", os.path.join(d, "o")])
            rank = open(os.path.join(d, "o", "pair_ranking.tsv")).read().splitlines()
            self.assertTrue(rank[1].startswith("Gallus gallus"))

    def test_analyse(self):
        r = analyse_protein(REF, [("x [A b]", REF)])
        self.assertEqual(r["rows"][0]["cmp"]["identity"], 1.0)


if __name__ == "__main__":
    unittest.main()
