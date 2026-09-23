# PROJECT 2 - OROPOUCHE STRUCTURE: LOCKED GATES
Locked: 2026-09-23 13:30 IST (Asia/Calcutta), before any outcome data is examined.
Builder 3 of 25. Repo: oropouche-structure (private, URL pending via parent).

## Definitions (locked)
- OUTBREAK set: OROV sequences sampled 2022-01-01 through 2026-12-31.
- HISTORICAL set: OROV sequences sampled before 2020-01-01. 2020-2021 excluded (buffer, reassortant emergence window 2010s-2021 per Nat Med 2024).
- Reference frames: L protein = RdRp ORF on L segment; N = nucleoprotein ORF on S; Gn/Gc = mature products of M polyprotein (cleavage per UniProt/literature boundaries).

## Success gates
- G1 (data): >=50 COMPLETE L-segment sequences byte-locked with SHA-256 checksums; M and S sets byte-locked alongside. Manifest = accession list + per-file checksums + fetch timestamps + source URLs.
- G2 (structure): structural models only from real sources (AlphaFold DB / experimentally grounded homology against PDB). Per-residue pLDDT tracked and reported; residues pLDDT<70 marked THIN on every figure and in every claim; no structural claim rests on thin regions. Positive control: model/alignment pipeline must recover the known Gc head domain architecture of experimental structure PDB 6H3X (TM-region excluded) and the catalytic endonuclease/RdRp motifs in L.
- G3 (mutation map): outbreak-vs-historical amino acid substitution list for L, N, Gn/Gc reproducible end-to-end from the manifest (re-run yields byte-identical tables). Reported with counts, frequencies, and segment lineage context (BR2015-2024 reassortant L/S vs M ancestry).
- G4 (validation before novel claims):
  PC1: motif scanner recovers canonical bunyavirus RdRp motifs (incl. motif C SDD) in 100% of complete L ORFs.
  PC2: glycosylation sequencer recovers the known N-glycan sequon present in PDB 6H3X from the reference M polyprotein.
  PC3: spike-in control - 5 engineered substitutions inserted into one historical L sequence are recovered exactly (5/5, no false positives) by the mutation pipeline.
  All three PCs must pass BEFORE any outbreak-vs-historical result is read.

## Negative-result protocol
If a gate fails, the failure is reported as-is in the paper. No re-fishing, no silent gate edits. Any gate change requires parent sign-off and is logged below.
## Gate change log
(none)
