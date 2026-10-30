# Recorded workflow and executable commands

Run from the package root with the selected Python environment active. Copy and edit `config/example.json` to `config/local.json` first. `run.py` sets configuration for each child process; it never uploads data or installs tools. Use a fresh run directory; preserve failed checkpoints for inspection. Some stages recycle only their own derived scratch files after checks. Do not point the run root at original inputs or a prior analysis directory.

## Reference-free candidate discovery and held-out association

```sh
python scripts/run.py --config config/local.json prepare --manifests
python scripts/run.py --config config/local.json discovery-counts
python scripts/run.py --config config/local.json discovery-screen
python scripts/run.py --config config/local.json heldout-counts
python scripts/run.py --config config/local.json heldout-tests
```

`discovery-counts` explicitly calls `run_qc_kmc.py --stage full` on the discovery manifest. The script's historical pilot default is not used. QC: cutadapt `-j 2 -a AAGTCGGAGGCCAAGCGGTC -A AAGTCGGATCGTAGCCATGT -O 8 -e 0.1 -q 20 --max-n 0.05 -m 50 --pair-filter=any`. Counting: KMC `-k31 -m2 -sm -t4 -ci2 -cs65535 -hp`; canonical counting is the default (no `-b`). Discovery uses presence >=2, >=6/8 in the focal sex and <=1/8 in the opposite sex. Frozen candidate list hashes are retained before held-out statistics. Candidate intersections preserve the right-hand individual count with `-ocright`.

The held-out test is two-sided Fisher exact, Bonferroni across all frozen F+M candidates (7,178,221 in the recorded analysis), corrected P <=0.05, with discovery-consistent direction. Exact binomial intervals are pointwise descriptive intervals. No association is adjusted for pedigree, missing batch or local depth. No reference file is read by these five stages.

## Reference placement, opportunity normalization and regional sensitivity

```sh
python scripts/run.py --config config/local.json prepare --reference
python scripts/run.py --config config/local.json exact-mapping
python scripts/run.py --config config/local.json genmap-toy
python scripts/run.py --config config/local.json reference-opportunities
python scripts/run.py --config config/local.json normalized-density
python scripts/run.py --config config/local.json existing-depth
python scripts/run.py --config config/local.json regional-sensitivity
python scripts/run.py --config config/local.json regional-sensitivity-figure
```

Exact mapping uses all reference sequences and both strands without mismatches. GenMap commands retained in the source are `index -A divsufsort --version-check 0` and `map -K 31 -E 0 -T 8 -fs -r --version-check 0`; its synthetic check precedes the full-reference run. Reference opportunity counts and candidate counts use start-based, 1-based inclusive windows. Regional panels are fixed before count thresholds 2/3/5/10 are compared; this is not repeated association testing or read-depth adjustment. `existing-depth` requires an existing same-cohort local depth matrix and its scope document; this code does not recover its upstream alignment/scaling commands.

## Figure 3 (final numbering)

```sh
python scripts/run.py --config config/local.json figure3-initial
python scripts/run.py --config config/local.json figure3-profiles
python scripts/run.py --config config/local.json figure3-discovery-panels
python scripts/run.py --config config/local.json figure3-all200
python scripts/run.py --config config/local.json figure3-check
```

These are the recorded display-preparation sequence, followed by the final all-200 renderer. Intermediate displays are not new manuscript figures. The fixed individual panels have 32,753 F-associated and 4,507 M-associated k-mers in Chr19:2,710,000–2,890,000. Filled versus open symbols distinguish held-out and discovery individuals. The all-200 display is descriptive; it is not a 200-individual validation cohort. Numeric assertions retain the recorded figure inputs. Individual source tables must be supplied/generated locally and are excluded from the candidate archive.

## Local genotype correlations, sequence predictions and secondary annotation

```sh
python scripts/run.py --config config/local.json regional-vcf
python scripts/run.py --config config/local.json local-correlations
python scripts/run.py --config config/local.json candidate-sequences
python scripts/run.py --config config/local.json matched-and-genotype-tests
python scripts/run.py --config config/local.json local-correlation-figure
python scripts/run.py --config config/local.json secondary-annotations
```

The regional correlation route masks DP <5/GQ <20 (DP >=3 sensitivity), requires per-sex call rates >=90% and MAF >=5%, excludes constant/undefined dosages, uses pairwise completeness >=90%, and bins pair summaries by 10 kb. No population LD decay model is fitted. The supplied matched comparison uses common sites/pairs in the entire 600-kb window. It is not the separate final SDR-restricted 193-site / 13,191-pair comparison; see the source-code gap in COVERAGE.md. The final 2x3 exact test applies Bonferroni only to QC-eligible coding/splice records in the two prespecified genes (28 records in the recorded run). Exploratory protein-guided NG86 estimates exclude incomplete codons and test 0/2/5-aa gap flanks. These are reference-segment comparisons, not branch-selection tests, phased Z/W proteins or functional proof.

`secondary-annotations` selects the largest unique-k-mer cluster for each direction on Chr34 and Contig6 and intersects gene bodies. It requires GFF and functional inputs. It does not rerun GWAS or test new candidate-gene significance.

## Synthetic checks only

```sh
python tests/smoke.py
python scripts/analyze_heldout_kmers.py --self-test
python scripts/map_heldout_kmers.py --self-test
```

The synthetic check creates only artificial inputs and derived outputs in the ignored `run/smoke/` area. It does not fetch study data, install tools, or execute external read counting. See `validation/test_report.json` for what was actually run for this candidate.
