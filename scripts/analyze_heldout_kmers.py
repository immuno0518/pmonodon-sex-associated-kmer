from release_config import *
"""Audit frozen candidate counts, aggregate held-out prevalence, and test association.

All original per-individual count databases are read-only inputs. Derived binary
databases count individuals, not reads. No individual information leaves this PC.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import os
import pathlib
import subprocess
import sys
import time

# Use isolated, pinned wheels rather than modifying the existing Anaconda setup.
import numpy as np
from scipy.stats import beta, fisher_exact, mannwhitneyu

ROOT = RUN_ROOT
OUT = ROOT / 'outputs/kmer'
BASE = ROOT / 'work/kmer/heldout_validation'
WORK = ROOT / 'work/kmer/heldout_statistics'
RESULT = OUT / 'heldout_statistics'
TOOL = KMC_TOOLS_BIN
STATUS = OUT / 'heldout_statistics_status.json'
SUFFIXES = ('.kmc_pre', '.kmc_suf')


def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = pathlib.Path(path)
    temp = path.with_suffix(path.suffix + '.partial')
    temp.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temp.replace(path)


def status(phase, **details):
    atomic_json(STATUS, dict(phase=phase, updated_at=time.strftime('%Y-%m-%d %H:%M:%S'), **details))


def run(args, log):
    with pathlib.Path(log).open('w', encoding='utf-8') as f:
        subprocess.run([str(TOOL), '-t2', '-hp'] + list(map(str, args)),
                       stdout=f, stderr=subprocess.STDOUT, check=True)


def db_hashes(prefix):
    return {suffix: digest(str(prefix) + suffix) for suffix in SUFFIXES}


def delete_derived_db(prefix):
    prefix = pathlib.Path(prefix).resolve()
    assert WORK.resolve() in prefix.parents, 'Refuse deletion outside statistics scratch space'
    for suffix in SUFFIXES:
        path = pathlib.Path(str(prefix) + suffix)
        if path.exists():
            path.unlink()


def proportion_ci(x, n):
    return (0.0 if x == 0 else float(beta.ppf(.025, x, n-x+1)),
            1.0 if x == n else float(beta.ppf(.975, x+1, n-x)))


def association(f, m, nf, nm, tests):
    odds, p = fisher_exact([[f, nf-f], [m, nm-m]], alternative='two-sided')
    return float(odds), float(p), min(1.0, float(p) * tests)


def audit(rows, frozen):
    audited = []
    qc_rows = []
    for i, row in enumerate(rows):
        sid = row['sample_id']
        d = BASE / sid
        done = json.loads((d / 'complete.json').read_text())
        assert done['sample_id'] == sid and done['sex'] == row['sex']
        assert done['frozen_candidates'] == frozen
        assert done['per_sample_presence_threshold'] == 2
        assert len(done['candidate_database_sha256']) == 4
        for name, expected in done['candidate_database_sha256'].items():
            assert name in [sex + '_candidate_counts' + suffix for sex in ('F', 'M') for suffix in SUFFIXES]
            assert digest(d / name) == expected, 'Checksum mismatch: ' + sid + '/' + name
        qc = json.loads((d / 'cutadapt.json').read_text())
        kmc = json.loads((d / 'kmc.json').read_text())['Stats']
        hist_totals = {}
        for group in ('F', 'M'):
            with (d / (group + '_candidate_count_histogram.tsv')).open() as f:
                hist = [tuple(map(int, line.split())) for line in f if line.strip()]
            assert all(count >= 2 and number >= 0 for count, number in hist)
            hist_totals[group] = sum(number for count, number in hist)
        qc_rows.append(dict(sample_id=sid, sex=row['sex'], sequencing_batch=row['sequencing_batch'],
                            output_read_pairs=qc['read_counts']['output'],
                            output_bases=qc['basepair_counts']['output'],
                            unique_counted_kmers=kmc['#Unique_counted_k-mers'],
                            total_kmer_observations=kmc['#Total no. of k-mers'],
                            F_candidate_present=hist_totals['F'], M_candidate_present=hist_totals['M']))
        audited.append(dict(sample_id=sid, sex=row['sex'], sha256=done['candidate_database_sha256']))
        if (i+1) % 10 == 0 or i+1 == len(rows):
            status('auditing_candidate_checksums', verified_samples=i+1, total=len(rows))
            print('Audited', i+1, '/', len(rows), flush=True)
    with (RESULT / 'sample_QC_and_candidate_burden.tsv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(qc_rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(qc_rows)
    identity = dict(frozen_candidates=frozen, samples=audited)
    atomic_json(RESULT / 'input_checksum_audit.json', identity)
    return identity, qc_rows


def aggregate(group, sample_sex, rows, identity_sha):
    selected = [r for r in rows if r['sex'] == sample_sex]
    d = WORK / (group + '_candidates_' + sample_sex + '_samples')
    d.mkdir(exist_ok=True)
    checkpoint = d / 'checkpoint.json'
    previous = None
    start = 0
    if checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        assert saved['input_identity_sha256'] == identity_sha
        start = saved['completed_samples']
        assert saved['sample_ids'] == [r['sample_id'] for r in selected[:start]]
        previous = d / saved['database']
        assert db_hashes(previous) == saved['database_sha256']
    for i in range(start, len(selected)):
        row = selected[i]
        status('aggregating_prevalence', candidate_group=group, sample_sex=sample_sex,
               completed_samples=i, total_samples=len(selected), current_sample=row['sample_id'])
        source = BASE / row['sample_id'] / (group + '_candidate_counts')
        target = d / ('acc_%03d' % i)
        binary = d / ('binary_%03d' % i) if previous else target
        # A candidate present with raw count >=2 contributes exactly one individual.
        run(['transform', source, '-ci2', 'set_counts', '1', binary], d / ('binary_%03d.log' % i))
        if previous:
            run(['simple', previous, binary, 'union', target, '-ocsum', '-cs255'], d / ('union_%03d.log' % i))
        atomic_json(checkpoint, dict(input_identity_sha256=identity_sha, completed_samples=i+1,
                    sample_ids=[r['sample_id'] for r in selected[:i+1]], database=target.name,
                    database_sha256=db_hashes(target)))
        if previous:
            delete_derived_db(previous); delete_derived_db(binary)
        previous = target
        if (i+1) % 10 == 0 or i+1 == len(selected):
            print(group, sample_sex, 'prevalence', i+1, '/', len(selected), flush=True)
    text = d / 'prevalence.sorted.tsv'
    if not text.exists():
        partial = text.with_suffix('.tsv.partial')
        run(['transform', previous, 'dump', '-s', partial], d / 'dump.log')
        partial.replace(text)
    return text


class SortedCounts:
    def __init__(self, path):
        self.handle = pathlib.Path(path).open()
        self.previous = ''
        self.current = self.next()

    def next(self):
        line = self.handle.readline()
        if not line:
            return None
        key, value = line.split()
        assert key > self.previous
        self.previous = key
        return key, int(value)

    def at(self, key):
        assert self.current is None or self.current[0] >= key, 'Unexpected k-mer outside frozen candidates'
        if self.current is not None and self.current[0] == key:
            value = self.current[1]
            self.current = self.next()
            return value
        return 0

    def finish(self):
        assert self.current is None
        self.handle.close()


def summarize(rows, frozen, texts, qc_rows):
    sizes = collections.Counter(r['sex'] for r in rows)
    nf, nm = sizes['F'], sizes['M']
    tests = sum(x['candidate_count'] for x in frozen.values())
    cache = {}
    patterns = collections.Counter()
    totals = {}
    all_path = RESULT / 'all_frozen_candidates_heldout_tests.tsv.gz'
    sig_path = RESULT / 'direction_consistent_bonferroni_associations.tsv.gz'
    fasta_path = RESULT / 'direction_consistent_associated_kmers.fa'
    header = 'candidate_group\tkmer\tdiscovery_focal_present_of_8\theldout_F_present\theldout_M_present\theldout_nF\theldout_nM\todds_ratio_F_vs_M\tfisher_two_sided_p\tbonferroni_p_all_candidates\tdirection_consistent\n'
    with gzip.open(str(all_path)+'.partial', 'wt', compresslevel=1) as all_out, \
         gzip.open(str(sig_path)+'.partial', 'wt', compresslevel=1) as sig_out, \
         pathlib.Path(str(fasta_path)+'.partial').open('w') as fasta:
        all_out.write(header); sig_out.write(header)
        for group in ('F', 'M'):
            status('testing_frozen_candidates', candidate_group=group, total_tests=tests)
            source = ROOT / 'work/kmer/discovery' / (group + '_associated_candidates')
            original_dump = WORK / (group + '_frozen.original.tsv')
            run(['transform', source, 'dump', original_dump], WORK / (group + '_frozen_verify.log'))
            assert digest(original_dump) == frozen[group]['sha256']
            sorted_dump = WORK / (group + '_frozen.sorted.tsv')
            run(['transform', source, 'dump', '-s', sorted_dump], WORK / (group + '_frozen_sorted.log'))
            fc, mc = SortedCounts(texts[group]['F']), SortedCounts(texts[group]['M'])
            n, sig, consistent_sig, perfect, reverse_sig = 0, 0, 0, 0, 0
            last = ''
            with sorted_dump.open() as candidates:
                for line in candidates:
                    key, discovery = line.split(); discovery = int(discovery)
                    assert key > last and len(key) == 31 and 6 <= discovery <= 8
                    last = key
                    f, m = fc.at(key), mc.at(key)
                    assert 0 <= f <= nf and 0 <= m <= nm
                    pair = (f, m)
                    if pair not in cache:
                        cache[pair] = association(f, m, nf, nm, tests)
                    odds, p, corrected = cache[pair]
                    consistent = (f/nf > m/nm) if group == 'F' else (m/nm > f/nf)
                    significant = corrected <= .05
                    out = '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%d\n' % (group,key,discovery,f,m,nf,nm,odds,p,corrected,consistent)
                    all_out.write(out)
                    if significant:
                        sig += 1
                        if consistent:
                            consistent_sig += 1; sig_out.write(out)
                            fasta.write('>%s_%07d F_present=%d M_present=%d\n%s\n' % (group,consistent_sig,f,m,key))
                        else:
                            reverse_sig += 1
                    is_perfect = (f == nf and m == 0) if group == 'F' else (m == nm and f == 0)
                    if is_perfect:
                        perfect += 1
                    patterns[(group,f,m)] += 1; n += 1
                    if n % 500000 == 0:
                        print(group, 'tested', n, flush=True)
            fc.finish(); mc.finish()
            assert n == frozen[group]['candidate_count']
            totals[group] = dict(tested=n, bonferroni_significant=sig,
                                 direction_consistent_bonferroni_significant=consistent_sig,
                                 reversed_direction_bonferroni_significant=reverse_sig,
                                 perfect_separation_in_heldout=perfect)
    for path in (all_path, sig_path, fasta_path):
        pathlib.Path(str(path)+'.partial').replace(path)
    with (RESULT/'contingency_pattern_summary.tsv').open('w') as out:
        out.write('group\tF_present\tM_present\tkmer_count\tfisher_p\tbonferroni_p\tfocal_sensitivity\tsensitivity_CI_low\tsensitivity_CI_high\tfocal_specificity\tspecificity_CI_low\tspecificity_CI_high\n')
        for (group,f,m), number in sorted(patterns.items()):
            focal, n_focal, other, n_other = (f,nf,m,nm) if group=='F' else (m,nm,f,nf)
            low_s, high_s = proportion_ci(focal,n_focal)
            low_c, high_c = proportion_ci(n_other-other,n_other)
            _, p, corrected = cache[(f,m)]
            out.write('%s\t%d\t%d\t%d\t%.12g\t%.12g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\n' % (group,f,m,number,p,corrected,focal/n_focal,low_s,high_s,(n_other-other)/n_other,low_c,high_c))
    qc = {}
    for metric in ('output_read_pairs','output_bases','unique_counted_kmers','F_candidate_present','M_candidate_present'):
        a = [r[metric] for r in qc_rows if r['sex']=='F']; b = [r[metric] for r in qc_rows if r['sex']=='M']
        qc[metric] = dict(F_median=float(np.median(a)),M_median=float(np.median(b)),
                          descriptive_mann_whitney_two_sided_p=float(mannwhitneyu(a,b,alternative='two-sided').pvalue))
    summary = dict(phase='sample_split_association_tests_complete_pending_genomic_interpretation',
                   nF=nf,nM=nm,tests=tests,alpha=.05,bonferroni_threshold=.05/tests,
                   candidate_groups=totals,QC_descriptive_comparison=qc,
                   method='Two-sided Fisher exact tests on held-out presence (raw count >=2), Bonferroni over ALL frozen F and M candidates; replication requires discovery-consistent direction.',
                   limitations=[
                       'Distinct held-out individuals, but a single full-sib family; not an independent population or family replication.',
                       'Nominal sample-level Fisher tests assume exchangeability; pedigree structure and unrecorded sequencing batches are not adjusted.',
                       'Bonferroni addresses multiplicity across correlated k-mers, not sample relatedness or confounding.',
                       'Absence means QC read count <2, not biological absence; raw candidate counts remain retained.',
                       'Perfect separation and confidence intervals are descriptive in this family; intervals are pointwise, not simultaneous.',
                       'Associated k-mers are not universal sex-specific markers, W-linked proof, or causal loci; mapping and genomic interpretation remain required.'],
                   sources=['https://docs.scipy.org/doc/scipy-1.5.2/reference/generated/scipy.stats.fisher_exact.html',
                            'https://github.com/refresh-bio/KMC'])
    atomic_json(RESULT/'association_summary.json',summary)
    atomic_json(RESULT/'result_checksums.json',{p.name:digest(p) for p in (all_path,sig_path,fasta_path)})
    return summary


def self_test():
    # Known probabilities, symmetry, all-absent handling, and exact binomial limits.
    assert abs(association(8,0,8,8,1)[1] - 2/math.comb(16,8)) < 1e-12
    assert association(0,0,8,8,10)[1:] == (1.0,1.0)
    assert abs(association(7,1,8,8,1)[1] - association(1,7,8,8,1)[1]) < 1e-12
    assert proportion_ci(0,92)[0] == 0 and proportion_ci(92,92)[1] == 1
    print('Statistical smoke checks passed.')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--self-test',action='store_true'); args=ap.parse_args()
    if args.self_test:
        self_test(); return
    WORK.mkdir(exist_ok=True); RESULT.mkdir(exist_ok=True)
    lock=WORK/'runner.lock'
    with lock.open('x') as f: f.write(str(os.getpid()))
    try:
        finished=json.loads((OUT/'heldout_validation_status.json').read_text())
        assert finished['phase']=='candidate_counts_complete_pending_statistical_analysis'
        assert len(finished['completed'])==184 and not (BASE/'runner.lock').exists()
        rows=list(csv.DictReader((OUT/'heldout_validation_manifest_184.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
        discovery=list(csv.DictReader((OUT/'discovery_manifest_16.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
        assert len(rows)==184 and len({r['sample_id'] for r in rows})==184
        assert not {r['sample_id'] for r in rows} & {r['sample_id'] for r in discovery}
        assert collections.Counter(r['sex'] for r in rows)=={'F':92,'M':92}
        frozen=json.loads((BASE/'frozen_candidates.json').read_text())
        for group in ('F','M'):
            assert digest(OUT/(group+'_candidate_kmers_discovery.tsv'))==frozen[group]['sha256']
        plan=dict(k=31,presence_threshold=2,frozen_candidates=frozen,heldout_nF=92,heldout_nM=92,
                  test='two-sided Fisher exact',multiplicity='Bonferroni over all frozen candidates',
                  alpha=.05,replication_direction='must agree with discovery sex',
                  no_universal_marker_or_W_assignment=True)
        plan_path=RESULT/'analysis_plan.json'
        if plan_path.exists(): assert json.loads(plan_path.read_text())==plan
        else: atomic_json(plan_path,plan)
        identity,qc_rows=audit(rows,frozen)
        identity_sha=digest(RESULT/'input_checksum_audit.json')
        texts={group:{sex:aggregate(group,sex,rows,identity_sha) for sex in ('F','M')} for group in ('F','M')}
        summary=summarize(rows,frozen,texts,qc_rows)
        status(summary['phase'], validation_completed=184, summary=str(RESULT/'association_summary.json'))
        print(json.dumps(summary),flush=True)
    except Exception as exc:
        status('failed',error=str(exc)); raise
    finally:
        lock.unlink()


if __name__=='__main__':
    main()
