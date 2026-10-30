# Scope and remaining reproducibility gaps

Main figures are numbered 1–7; supplementary figures S1–S12. The source-to-revision crosswalk below follows the final k-mer-main revision. It is not a statement that every original analysis is included.

| Final item / analysis | Candidate contents | Limitation |
|---|---|---|
| Figure 3 A–C | Recorded mapping/normalization and display scripts; three aggregate panel TSV files | Original reference, candidate databases and frequency array excluded; original placement not rerun |
| Figure 3 D/E | Recorded discovery-count extraction, all-200 display and readback scripts | Individual source tables/sex mappings excluded; cannot rerender these panels from included tables alone |
| Figure S7 | Regional VCF access, local dosage correlations, full-window common-pair comparison and plot scripts | Requires private 200-individual VCF and explicit sample labels; the separate final SDR-restricted matched comparison (193 sites / 13,191 pairs) has summary evidence, but its generating script was not recovered in the audited directories |
| Figure S8 | Fixed-panel threshold sensitivity and figure scripts | Requires retained individual candidate databases and existing local-depth inputs |
| Figure S9 | Exact reference-unique opportunity calculation and normalization | Historical supplementary renderer (`finalize_kmer_regional_assessment.py`) not included; output calculation is covered, exact source layout not recreated |
| Tables S4/S5 | Bounded secondary gene-body intersection and annotation lookup | Requires corresponding GFF/functional annotations; no original annotation data distributed |
| Tables S6/S7 | Recorded CDS consequence, exact 2x3 tests and exploratory NG86 code | Reference, prior CDS and phased-table prerequisites excluded; full candidate-sequence stage not rerun |
| Figure S10 | Final common-rule vendor gene-body comparison | Source/result audit exists elsewhere; its analysis/plot script is not recovered from this release's three audited source directories |
| Figures 1/2/4/5/6/7; S1–S6/S11/S12 | See final revision for original sources | Original GWAS/ED/FST/segregation, phylogenetic/domain, RNA, PCR/Sanger, KASP and legacy-figure pipelines are not represented as recovered here |
| Short-read candidate allele audits | Methods described in revision | BAM/pileup/long-read competitive or structural-window audit scripts not included in this candidate; no full structural breakpoint or phased W sequence claim |
| RNA quantification | Original report says HISAT2 2.2.1 and RSEM 1.3.3, reported defaults | Per-run invocation logs not recovered; this candidate does not provide a reconstructed RNA alignment pipeline |
| TE comparison | Excluded | Not verified complete; no enrichment/absence conclusion |

Previous main figures 3/4/5/6 became main figures 4/5/6/7. Previous supplementary S1/S2/S3/S4/S5/S6 became S2/S3/S4/S5/S6/S11. Previous revision S10/S11/S12/S13 became S8/S9/S7/S10; previous R1 S8 became S12.

Packaged tests validate defined functions, small synthetic data flows and integrity only. They cannot establish reproduction of millions of original candidate tests, the original reference's unique-start denominator or any biological validation. Frozen-candidate association, reference placement, descriptive correlations, sequence predictions and functional validation remain separate evidence levels.

`finalize_reviewer_statistics.py` retains the source's common-site comparison across the entire extracted 600-kb window. Its output `matched_pair_group_comparison.json` must not be relabelled as the final SDR-restricted `matched_SDR_pair_comparison.json`. Recover the latter's original generating code before claiming that the 193-site / 13,191-pair calculation is provided here.

Before using this package as the complete response to the public-code request, confirm whether the uncovered custom analyses also need recovered source code. Do not state that all custom manuscript scripts or all figures are publicly reproducible from this candidate alone.
