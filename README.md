# ProteinSequenceAlignment
For genes that have not been resequenced, might have missing nucleotides in exons or truncated terminals

## protaln: bird TF / partner protein mining, alignment, QC and scoring

Pure-Python (3.8+, no dependencies) pipeline that, for a transcription factor and its interacting partner:

1. **Mines** bird (taxid 8782) protein sequences from NCBI via E-utilities (`--tf-gene`, `--partner-gene`), or reads local FASTA (`--tf-fasta`, `--partner-fasta`) e.g. Vertebrate Genomes Project translations.
2. **Aligns** each species to the known reference (BLOSUM62, free end gaps) and reports per-residue substitutions.
3. **Detects motifs** (NLS, zinc finger, LxxLL, bHLH, N-glycosylation by default).
4. **QC**: length, internal stops, ambiguous residues, missing start Met.
5. **Consensus** across QC-passing species, with disagreement regions vs the reference/variable columns.
6. **Completeness**: reference coverage, N/C-terminal truncation, internal gaps.
7. **Scores** per species and a `pair_ranking.tsv` (pair score = the weaker of TF/partner) to choose a construct for the transfection plasmid.

```
python -m protaln --tf-ref tf_ref.fa --partner-ref partner_ref.fa \
    --tf-gene SOX2 --partner-gene POU5F1 --outdir results
python -m unittest discover -s tests
```
Outputs in `--outdir`: `{tf,partner}_sequences.fa`, `*_per_species.tsv`, `*_consensus.tsv`, `*_disagreement_regions.tsv`, `pair_ranking.tsv`.
