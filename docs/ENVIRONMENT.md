# Software records and testing environment

| Component | Recorded study/revision evidence | Status |
|---|---|---|
| Python | 3.8.5, from held-out `analysis_environment_and_checks.json` | Recorded association environment |
| KMC | 3.2.4, from `software_versions.json` and command strings | Recorded counting version |
| cutadapt | 4.4, from the same record and installed metadata | Recorded QC version |
| NumPy / SciPy | 1.24.4 / 1.10.1, from held-out environment record | Recorded association versions |
| pyahocorasick | 2.1.0, retained sequence-mapping dependency metadata | Recovered installed version; no separate per-run version invocation |
| Biopython | 1.83, retained regional-analysis dependency metadata | Recovered installed version; no separate per-run version invocation |
| GenMap | 1.3.0, specified in the reference-opportunity script | Recorded configured version; no new full-reference run in this release |
| dnaio / isal / xopen / zlib-ng | 1.2.2 / 1.7.2 / 2.0.2 / 0.5.1 | Recorded QC dependency list |
| psutil / matplotlib | Exact historical versions not recovered | Required by resource-monitoring / plotting scripts; do not infer historical versions from release tests |

The original analysis used separate dependency directories; a single complete historical lockfile was not recovered. `requirements-recorded.txt` lists recovered pins and identifies the remaining unpinned tools. It is not a claim that a clean historical environment has been reconstructed. Supply KMC 3.2.4 and GenMap 1.3.0 separately, subject to their licences. No binaries or third-party libraries are bundled.

The candidate's actual synthetic-test environment and results are saved in `validation/test_report.json`. The Python 3.12 test environment is distinct from the original Python 3.8.5 environment. It establishes only the checks listed in that report. Use a separate environment for the old dependencies if reproducing the recorded workflow; do not replace an existing analysis environment to run these tests.

For current smoke-test dependencies, `requirements-smoke.txt` records the versions used for this candidate. No cutadapt, KMC or GenMap computation is claimed by those Python-only tests. No full-cohort, full-reference, RNA alignment, SNP calling or TE calculation was performed during candidate packaging.

Upstream software: [KMC](https://github.com/refresh-bio/KMC), [cutadapt](https://github.com/marcelm/cutadapt), [GenMap](https://github.com/cpockrandt/genmap), [SciPy](https://scipy.org/), [Biopython](https://biopython.org/). These are software links, not a public address for this study's code.
