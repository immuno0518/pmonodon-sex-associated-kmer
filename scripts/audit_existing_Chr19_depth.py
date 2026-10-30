from release_config import *
"""Audit pre-existing same-cohort depth; never treat it as independent validation."""
import csv
import hashlib
import json
import statistics
import sys
from collections import defaultdict, Counter
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
import numpy as np
from scipy.stats import mannwhitneyu,spearmanr


def main():
    source=input_path('depth_matrix')
    scope=input_path('depth_scope')
    with source.open(encoding='utf-8-sig') as f:
        depth=list(csv.DictReader(f,delimiter='\t'))
    with (ROOT/'outputs/kmer/heldout_statistics/sample_QC_and_candidate_burden.tsv').open() as f:
        heldout=list(csv.DictReader(f,delimiter='\t'))
    grouped=defaultdict(list)
    for r in depth:
        grouped[r['sample_id']].append(r)
    assert len(depth)==12000 and len(grouped)==200
    assert Counter(r[0]['sex'] for r in grouped.values())=={'F':100,'M':100}
    for rows in grouped.values():
        assert len(rows)==60 and len({r['sex'] for r in rows})==1
        assert len({(r['chrom'],r['start'],r['end']) for r in rows})==60
        assert abs(statistics.median(float(r['mean_depth']) for r in rows)-float(rows[0]['sample_median_depth']))<0.00011
        assert all(abs(float(r['mean_depth'])/float(r['sample_median_depth'])-float(r['normalized_depth']))<0.000001 for r in rows)
    joined=[]
    for sample in heldout:
        sid=sample['sample_id']
        rows=grouped[sid]
        assert rows[0]['sex']==sample['sex']
        selected=[]
        for r in rows:
            weight=max(0,min(int(r['end']),2890000)-max(int(r['start']),2710000)+1)
            if weight:
                selected.append((r,weight))
        assert sum(w for r,w in selected)==180001
        r=dict(sample_id=sid,sex=sample['sex'],target_mean_depth_10kb_approx=sum(float(r['mean_depth'])*w for r,w in selected)/180001,
               target_depth_relative_to_local_600kb_median=sum(float(r['normalized_depth'])*w for r,w in selected)/180001,
               median_depth_in_local_600kb=float(rows[0]['sample_median_depth']),
               output_read_pairs=int(sample['output_read_pairs']),F_frozen_candidate_present=int(sample['F_candidate_present']),M_frozen_candidate_present=int(sample['M_candidate_present']))
        joined.append(r)
    stats={}
    correlations=[]
    for measure in ('target_mean_depth_10kb_approx','target_depth_relative_to_local_600kb_median','median_depth_in_local_600kb'):
        groups={sex:[r[measure] for r in joined if r['sex']==sex] for sex in 'FM'}
        stats[measure]={sex:{'n':len(v),'q1':float(np.quantile(v,.25)),'median':float(np.median(v)),'q3':float(np.quantile(v,.75))} for sex,v in groups.items()}
        stats[measure]['descriptive_unadjusted_mann_whitney_p']=float(mannwhitneyu(groups['F'],groups['M'],alternative='two-sided').pvalue)
        for sex in 'FM':
            rows=[r for r in joined if r['sex']==sex]
            for burden in ('F_frozen_candidate_present','M_frozen_candidate_present','output_read_pairs'):
                rho,p=spearmanr([r[measure] for r in rows],[r[burden] for r in rows])
                correlations.append(dict(sex=sex,depth_measure=measure,whole_pool_burden_or_QC=burden,descriptive_rho=float(rho),descriptive_unadjusted_p=float(p)))
    out=ROOT/'outputs/kmer/reference_mappability'
    with (out/'existing_same_cohort_depth_heldout_join.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(joined[0]),delimiter='\t');w.writeheader();w.writerows(joined)
    audit=dict(created_at=datetime.now().astimezone().isoformat(),source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
               source_scope_document=str(scope),source_scope_sha256=hashlib.sha256(scope.read_bytes()).hexdigest(),
               existing_depth_not_newly_recomputed=True,independent_cohort=False,
               audit={'rows':12000,'individuals':200,'bins_per_individual':60,'heldout_joined':184,'heldout_sex_labels_consistent':True,'normalizer':'Median of the SAME 60 local 10kb bins spanning Chr19:2500001-3100000; not a genome-wide normalizer','ratio_recalculation':'passed'},
               target='Chr19:2710000-2890000',target_weighting='Overlap-weighted 10kb bin mean; not base-resolution depth. One boundary bin contributes 1bp.',
               descriptive_depth_comparison=stats,within_sex_correlations=correlations,
               limitations=['Source is pre-existing same-cohort evidence; provenance scope document does not establish independent validation or a newly rerun alignment-depth analysis.',
                            'Local-median normalization can absorb broad regional biological or technical differences; it is not a genome-wide sequencing-depth correction.',
                            'Burden correlations use the WHOLE frozen candidate pools, not a region-specific per-individual kmer panel.',
                            'These summaries do not adjust original kmer Fisher tests, prove equal local depth, or distinguish biological copy variation from alignment bias.',
                            'Comparable individual depth for Chr34 and Contig6 is not provided by this Chr19-only matrix.'])
    (out/'existing_same_cohort_depth_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'audit':audit['audit'],'depth_comparison':stats,'strongest_within_sex_correlations':sorted(correlations,key=lambda r:abs(r['descriptive_rho']),reverse=True)[:4]},indent=2))


if __name__=='__main__':
    main()
