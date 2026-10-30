# Changes from the recorded scripts

Each released analysis script has an original-file SHA-256 and a released-file SHA-256 in `provenance/script_sources.json`. Source IDs are stable aliases for the retained local source inventory, not claims of a public source location. Original files were left unchanged. No development history or original commit dates have been reconstructed.

| Files | Changes | Calculation retained |
|---|---|---|
| All 22 sourced scripts | Added shared configuration import, replaced machine-specific roots and private working-directory labels, removed injected local library search paths | Existing formulas, thresholds, sample-size gates and study result assertions |
| `run_qc_kmc.py`, `validate_heldout_kmers.py` | Use installed tools/current Python environment; UTF-8 input lists replace a Windows-specific encoding; no bundled dependency path | Adapter sequences; cutadapt and KMC settings; count preservation and scratch cleanup |
| KMC callers | Executables from configuration/PATH; two Windows executable checksum assertions become optional configured checksum checks | KMC operation order and counter modes |
| `map_heldout_kmers.py` | Reference path and reference dimensions from configuration; reference identity comes from a newly generated same-name sequence manifest | Exact both-strand matching, 1-based coordinates, unique/multiple/unmapped categories, per-sequence checksums |
| `server_reference_mappability.py`, `test_reference_genmap.py` | Use configured local reference and tool; keep input sequence IDs; write directly to the analysis input folder; condition CPU-affinity call on OS support; optional binary checksum replaces a stored archive check | GenMap k31/e0 and both-strand options, byte-frequency interpretation, synthetic validation and checkpoint checks |
| `analyze_reference_mappability.py`, `create_kmer_manuscript_figure.py`, `revise_kmer_manuscript_figure.py` | Replace private renamed-to-original chromosome key with an identity manifest, generated from the exact input FASTA | Reference start denominators, window sums and density arithmetic |
| `audit_existing_Chr19_depth.py` | Configure depth matrix and its provenance-scope document | Original overlap weighting and same-600-kb median normalization |
| `reviewer_genomics_common.py`, `analyze_reviewer_genes.py` | Configure reference, GFF, VCF source directory, previous CDS directory, phased-variant table and focal GFF IDs | Indexed extraction, genotype masks, reference reconstruction, individual-variant consequences and exploratory NG86 calculation |
| `secondary_gene_intersections.py` | Extracted only the bounded annotation section of the R2 integration audit; configured GFF, functional table and focal IDs; small file-I/O helpers replace manuscript-audit dependencies | Largest cluster by direction on each secondary sequence, gene-body overlap and annotation lookup |

The newly added `release_config.py`, `prepare_run.py`, `run.py`, synthetic checks and documentation are release utilities, not recovered original analysis code. They must not be described as original study scripts. The original scripts' module text and scientific limitations are retained except where a path or configuration reference required adjustment. Study-specific assertions make these scripts a documented study workflow rather than a general-purpose pipeline.

The private chromosome key is neither packaged nor recreated. All reference-dependent stages now require one consistent sequence-ID system in FASTA, GFF, VCF and derived tables. Identity-manifest generation records only IDs, lengths and sequence hashes from a user-supplied reference.

Two historical details are retained transparently: the original Fisher script cites an older SciPy documentation URL, while the recorded execution used SciPy 1.10.1; the local genotype function masks DP/GQ when those FORMAT fields are supplied. The release dispatcher requires both fields for the regional VCF preparation route. Source records, rather than a newly installed environment, define the historical software evidence.
