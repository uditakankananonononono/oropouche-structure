# oropouche/ - Project 2: Oropouche virus structural + mutational atlas (builder 3/25)

Byte-locked, validation-gated outbreak-vs-historical substitution catalog for OROV L, N, Gn/Gc with experimental/homolog-anchored structural context.

- `GATES_LOCKED.md` - gates + positive controls locked before outcome data (standalone first commit)
- `code/` - fetch_data.py, extract_cds.py, mutation_pipeline.py (PC3 built in, refuses to emit on failure), make_consensus.py
- `data/` - manifest.tsv (791 rows, per-record SHA-256), BYTELOCK.txt (set checksums), structure_targets.fasta. Bulk nucleotide FASTAs re-derivable from manifest accessions via code/fetch_data.py; checksums verify identity.
- `results/` - substitutions.tsv (36 rows, byte-identical on re-run), validation CSVs, structures/ (ESMFold PDBs - thin, claims restricted per gate G2)
- `figures/` - F1-F6 PNGs
- `paper/` - paper.tex + paper.pdf (17 pp)

Headline: Gc I1203T near-fixed in outbreak (97.3%); L polymerase-core cluster (787-856); N M170T 8.1 A from RNA in SBV 4JNG anchor. ESMFold single-sequence folding failed the 6H3X positive control (17.7 A global RMSD) - preserved negative result.
