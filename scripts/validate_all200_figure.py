from release_config import *
"""Independent readback of the all-200 plotting table and extracted discovery counts."""
import csv,gzip,json,hashlib
from pathlib import Path
from collections import Counter
from statistics import median
BASE=RUN_ROOT/'outputs/kmer'
OUT=BASE/'manuscript_main_figure_v3'
def rows(p):
    with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
data=rows(OUT/'panels_DE_all200_detection.tsv')
original={(r['sample_id'],r['candidate_group']):r for r in rows(BASE/'regional_sensitivity/per_individual_region_threshold_sensitivity.tsv')
          if r['region']=='Chr19_target' and r['read_count_threshold']=='2'}
discovery={r['sample_id']:r['sex'] for r in rows(BASE/'discovery_manifest_16.tsv')}
assert len(discovery)==16 and len({r['sample_id'] for r in data})==200
panels={r['kmer']:r['candidate_group'] for r in rows(BASE/'regional_sensitivity/fixed_regional_panels.tsv') if r['region']=='Chr19_target'}
counts={};nonzero_rows=0
for sid,sex in discovery.items():
    counts[sid]=Counter();seen=set()
    with gzip.open(OUT/(sid+'.discovery_Chr19_raw_counts.tsv.gz'),'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            assert r['kmer'] not in seen and panels[r['kmer']]==r['candidate_group']
            assert int(r['retained_read_count'])>=2
            counts[sid][r['candidate_group']]+=1;seen.add(r['kmer']);nonzero_rows+=1
for r in data:
    key=r['sample_id'],r['candidate_group']
    if r['cohort']=='Heldout':
        old=original[key]
        assert r['sex']==old['sex'] and int(r['detected_kmers'])==int(old['detected_kmers']) and int(r['panel_kmers'])==int(old['panel_kmers'])
    else:
        assert r['cohort']=='Discovery' and discovery[r['sample_id']]==r['sex']
        assert int(r['detected_kmers'])==counts[r['sample_id']][r['candidate_group']]
    assert abs(float(r['detected_percent'])-100*int(r['detected_kmers'])/int(r['panel_kmers']))<1e-10
summary=[]
for g in 'FM':
    for sex in 'FM':
        subset=[r for r in data if r['candidate_group']==g and r['sex']==sex]
        assert len(subset)==100 and len({r['sample_id'] for r in subset})==100
        assert Counter(r['cohort'] for r in subset)=={'Heldout':92,'Discovery':8}
        summary.append({'candidate_group':g,'sex':sex,'n':100,'discovery':8,'heldout':92,
                        'median_detection_percent':median(float(r['detected_percent']) for r in subset)})
for name in ['panel_A_chromosome_density.tsv','panel_B_whole_sequence_profiles.tsv','panel_C_Chr19_10kb_density.tsv']:
    assert sha(OUT/name)==sha(BASE/'manuscript_main_figure_v2'/name)
p=json.loads((OUT/'figure_provenance.json').read_text())
for name,h in p['figure_and_source_data_sha256'].items():assert sha(OUT/name)==h
report={'status':'passed','all200_table_rows':len(data),'discovery_nonzero_count_rows_reread':nonzero_rows,
        'heldout_original_values_unchanged':True,'ABC_source_data_unchanged':True,'output_hashes_verified':True,'group_summaries':summary}
(OUT/'all200_figure_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
