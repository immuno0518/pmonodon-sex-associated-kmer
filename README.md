# Penaeus monodon sex-associated k-mer and regional analyses

Code accompanying manuscript AQUACULTURE-D-26-04327. The rights holder confirmed ownership and authorized public availability on 2026-10-31. No open-source licence is granted. See [NOTICE.md](NOTICE.md).

Repository: [immuno0518/pmonodon-sex-associated-kmer](https://github.com/immuno0518/pmonodon-sex-associated-kmer).

The scripts are portable copies or bounded extractions of recorded project scripts. They cover read QC and k31 counting, discovery in 8 females and 8 males, frozen-candidate testing in the remaining 92 females and 92 males, exact reference placement, reference-unique start normalization, regional sensitivity, local genotype correlations and candidate-gene sequence summaries. All 200 resequenced individuals belong to one full-sib family. The held-out set is not an independent population or pedigree.

Discovery precedes reference placement and requires paired reads, not a reference assembly. Placement, reference-unique denominators and gene annotation require the corresponding reference or annotation files. The unpublished reference, GFF, reads, private chromosome-name key, individual sex-mapping tables and credentials are excluded. No TE analysis is distributed.

## Start

1. Install the tools and libraries in [ENVIRONMENT.md](docs/ENVIRONMENT.md).
2. Copy `config/example.json` to `config/local.json`; configure your authorized inputs and tool paths using [INPUTS.md](docs/INPUTS.md).
3. Run `python scripts/run.py --config config/local.json prepare` to create a separate run directory. Follow [COMMANDS.md](docs/COMMANDS.md) for the actual stage order.
4. Run `python tests/smoke.py` for synthetic checks without study data. These checks do not require KMC or GenMap. The separate recorded GenMap synthetic check requires GenMap.

Read [COVERAGE.md](docs/COVERAGE.md) before interpreting reproducibility. Full read counting, the original reference placement and all manuscript analyses have not been rerun for this package. Study-specific result assertions remain in several scripts and intentionally stop runs using incompatible data. Only aggregate Figure 3 A–C source tables are included; individual-level D/E inputs are excluded.

Script hashes and source identities appear in `provenance/script_sources.json`; [ADAPTATIONS.md](docs/ADAPTATIONS.md) explains portability changes. `MANIFEST.sha256` covers the delivered files. See [NOTICE.md](NOTICE.md) for rights status. The repository documents the supplied code and its coverage limits; it does not claim that every manuscript analysis or figure can be reproduced from the included files.
