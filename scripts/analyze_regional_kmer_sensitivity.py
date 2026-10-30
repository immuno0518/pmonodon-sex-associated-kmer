from release_config import *
"""Local descriptive sensitivity of fixed mapped panels using retained raw counts."""
import collections
import csv
import gzip
import hashlib
import json
import os
import statistics
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
import numpy as np
from scipy.stats import spearmanr

WORK=ROOT/'work/kmer/regional_sensitivity'
OUT=ROOT/'outputs/kmer/regional_sensitivity'
INPUT=ROOT/'work/kmer/heldout_validation'
MAPPING=ROOT/'outputs/kmer/heldout_mapping/all_replicated_candidates_reference_mapping.tsv.gz'
KMC=KMC_BIN
TOOLS=KMC_TOOLS_BIN
STATUS=ROOT/'outputs/kmer/regional_sensitivity_status.json'
THRESHOLDS=(2,3,5,10)
REGIONS={'Chr19_target':('Chr19',2710000,2890000),'Chr34':('Chr34',1,38034696),'Contig6':('Contig6',1,375397)}
SUFFIXES=('.kmc_pre','.kmc_suf')


def now():return datetime.now().astimezone().isoformat()


def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def save(p,v):
    temp=p.with_name(p.name+'.partial')
    temp.write_text(json.dumps(v,indent=2)+'\n',encoding='utf8');temp.replace(p)


def state(phase,**kw):
    save(STATUS,dict(phase=phase,pid=os.getpid(),updated_at=now(),**kw))


def tool(args,log):
    with log.open('w',encoding='utf8') as f:
        subprocess.run(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,check=True)


def db_hashes(p):return {s:digest(str(p)+s) for s in SUFFIXES}


def remove_scratch(p):
    p=p.resolve();assert WORK.resolve() in p.parents
    if p.exists():p.unlink()


def prepare_panels():
    meta=WORK/'panels.complete.json'
    panel_file=OUT/'fixed_regional_panels.tsv'
    if meta.exists():
        m=json.loads(meta.read_text())
        assert digest(MAPPING)==m['mapping_sha256'] and digest(panel_file)==m['panel_file_sha256']
        for g in 'FM':assert db_hashes(WORK/('panel_'+g))==m['panel_database_sha256'][g]
        with panel_file.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        return {r['kmer']:r for r in rows},m
    rows=[]
    with gzip.open(MAPPING,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if int(r['exact_reference_hits'])!=1:continue
            c,a=r['unique_chr'],int(r['unique_start_1based'])
            for name,(chrom,start,end) in REGIONS.items():
                if c==chrom and start<=a<=end:
                    rows.append(dict(region=name,candidate_group=r['candidate_group'],kmer=r['kmer'],heldout_F_present=r['heldout_F_present'],heldout_M_present=r['heldout_M_present'],perfect_heldout_separation=r['perfect_heldout_separation']))
    assert len(rows)==50915 and len({r['kmer'] for r in rows})==50915
    with panel_file.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    panel_hashes={}
    for g in 'FM':
        seqs={r['kmer'] for r in rows if r['candidate_group']==g}
        fa=WORK/('panel_'+g+'.fa')
        fa.write_text(''.join('>q'+str(i)+'\n'+s+'\n' for i,s in enumerate(sorted(seqs))))
        db=WORK/('panel_'+g)
        temp=WORK/('temp_'+g);temp.mkdir(exist_ok=True)
        tool([KMC,'-k31','-t2','-m2','-sm','-ci1','-cs65535','-hp','-fm',fa,db,temp],WORK/('panel_'+g+'.kmc.log'))
        dump=WORK/('panel_'+g+'.dump.tsv')
        tool([TOOLS,'-t2','-hp','transform',db,'-ci1','dump',dump],WORK/('panel_'+g+'.dump.log'))
        observed={}
        with dump.open() as f:
            for line in f:
                k,n=line.split();observed[k]=int(n)
        assert set(observed)==seqs and set(observed.values())=={1}
        panel_hashes[g]=db_hashes(db)
    expected={}
    for region in REGIONS:
        for g in 'FM':
            subset=[r for r in rows if r['region']==region and r['candidate_group']==g]
            expected[region+'_'+g]={'panel_kmers':len(subset),'F_prevalence_sum':sum(int(r['heldout_F_present']) for r in subset),'M_prevalence_sum':sum(int(r['heldout_M_present']) for r in subset)}
    m=dict(created_at=now(),mapping_sha256=digest(MAPPING),panel_file_sha256=digest(panel_file),panel_database_sha256=panel_hashes,expected=expected)
    save(meta,m)
    return {r['kmer']:r for r in rows},m


def check_input(sid):
    d=INPUT/sid;m=json.loads((d/'complete.json').read_text())
    for name,h in m['candidate_database_sha256'].items():assert digest(d/name)==h
    return m['candidate_database_sha256']


def sample(s,panels,metadata):
    sid=s['sample_id'];assert sid.replace('-','').isalnum()
    d=WORK/sid;d.mkdir(exist_ok=True)
    checkpoint=d/'complete.json';dest=OUT/(sid+'.regional_raw_counts.tsv.gz')
    if checkpoint.exists():
        m=json.loads(checkpoint.read_text())
        assert m['panel_file_sha256']==metadata['panel_file_sha256'] and digest(dest)==m['raw_panel_counts_sha256']
        return m
    before=check_input(sid)
    counters={(region,g):collections.Counter() for region in REGIONS for g in 'FM'}
    temp=dest.with_name(dest.name+'.partial')
    with gzip.open(temp,'wt',encoding='ascii') as f:
        f.write('region\tcandidate_group\tkmer\tretained_read_count\n')
        for g in 'FM':
            intersect=d/('intersection_'+g);dump=d/('intersection_'+g+'.tsv')
            # Only incomplete derived scratch can be recreated; original DBs are read-only.
            for suffix in SUFFIXES:remove_scratch(Path(str(intersect)+suffix))
            remove_scratch(dump)
            tool([TOOLS,'-t2','-hp','simple',WORK/('panel_'+g),'-ci1',INPUT/sid/(g+'_candidate_counts'),'-ci2','intersect',intersect,'-ocright','-cs65535'],d/(g+'.intersection.log'))
            tool([TOOLS,'-t2','-hp','transform',intersect,'dump',dump],d/(g+'.dump.log'))
            seen=set()
            with dump.open() as source:
                for line in source:
                    k,n=line.split();n=int(n)
                    assert k not in seen;seen.add(k)
                    r=panels[k];assert r['candidate_group']==g and 2<=n<=65535
                    counters[(r['region'],g)][n]+=1
                    f.write(r['region']+'\t'+g+'\t'+k+'\t'+str(n)+'\n')
    assert check_input(sid)==before
    temp.replace(dest)
    rows=[]
    for (region,g),hist in counters.items():
        total=metadata['expected'][region+'_'+g]['panel_kmers']
        for t in THRESHOLDS:
            present=sum(v for n,v in hist.items() if n>=t)
            rows.append(dict(sample_id=sid,sex=s['sex'],region=region,candidate_group=g,read_count_threshold=t,panel_kmers=total,detected_kmers=present,detected_fraction=present/total))
    m=dict(sample_id=sid,sex=s['sex'],completed_at=now(),panel_file_sha256=metadata['panel_file_sha256'],raw_panel_counts_sha256=digest(dest),input_database_sha256_before_and_after=before,source_files_unchanged=True,summary_rows=rows,
           count_histograms={region+'_'+g:dict(hist) for (region,g),hist in counters.items()})
    save(checkpoint,m)
    # Durable compressed raw panel counts and hashes replace only derived temporary intersections.
    for g in 'FM':
        for suffix in SUFFIXES:remove_scratch(Path(str(d/('intersection_'+g))+suffix))
        remove_scratch(d/('intersection_'+g+'.tsv'))
    return m


def summarize(results,metadata):
    rows=[r for m in results for r in m['summary_rows']]
    with (OUT/'per_individual_region_threshold_sensitivity.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    for name,expected in metadata['expected'].items():
        region,g=name.rsplit('_',1)
        for sex in 'FM':
            n=sum(r['detected_kmers'] for r in rows if r['region']==region and r['candidate_group']==g and r['sex']==sex and r['read_count_threshold']==2)
            assert n==expected[sex+'_prevalence_sum'],(name,sex,n,expected)
    depth={}
    with (ROOT/'outputs/kmer/reference_mappability/existing_same_cohort_depth_heldout_join.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):depth[r['sample_id']]=r
    summaries=[];correlations=[]
    for region in REGIONS:
        for g in 'FM':
            for threshold in THRESHOLDS:
                item=dict(region=region,candidate_group=g,read_count_threshold=threshold,panel_kmers=metadata['expected'][region+'_'+g]['panel_kmers'])
                for sex in 'FM':
                    selected=[r for r in rows if r['region']==region and r['candidate_group']==g and r['read_count_threshold']==threshold and r['sex']==sex]
                    values=[r['detected_fraction'] for r in selected];assert len(values)==92
                    item[sex]={'n':92,'min':min(values),'q1':float(np.quantile(values,.25)),'median':float(np.median(values)),'q3':float(np.quantile(values,.75)),'max':max(values)}
                    metrics=['output_read_pairs']
                    if region=='Chr19_target':metrics+=['target_mean_depth_10kb_approx','target_depth_relative_to_local_600kb_median']
                    for metric in metrics:
                        v=[float(depth[r['sample_id']][metric]) for r in selected]
                        rho=float(spearmanr(values,v).correlation) if len(set(values))>1 and len(set(v))>1 else None
                        correlations.append(dict(region=region,candidate_group=g,threshold=threshold,sex=sex,covariate=metric,n=92,descriptive_spearman_rho=rho))
                summaries.append(item)
    save(OUT/'regional_threshold_summary.json',dict(completed_at=now(),summary=summaries,within_sex_descriptive_correlations=correlations,
         validation='All regional count>=2 sex-wise sums exactly match previously mapped heldout per-kmer prevalences; all original source DB SHA256 unchanged',
         limitations=['Panels are conditional on unique reference matching and previously selected association; this is sensitivity, not a new independent association test.',
                      'Counts below 2 were not retained in original count databases; absence is censored count<2, not established biological absence.',
                      'Changing read-count threshold does not equalize sequencing depth.',
                      'Source counts saturate at 65535; raw count tables retain that value without interpreting it as uncapped.',
                      'Chr19 depth is pre-existing same-cohort 10kb approximation, locally normalized; comparable Chr34/Contig6 depth is unavailable.',
                      'Single full-sib family, missing batch covariate, and possible copy/haplotype/alignment effects remain; no adjusted or causal claim.']))
    save(OUT/'independent_prevalence_and_input_validation.json',dict(completed_at=now(),individuals=184,input_database_files_verified_before_and_after=736,all_sources_unchanged=True,panel_sizes={k:v['panel_kmers'] for k,v in metadata['expected'].items()},all_region_group_sex_prevalence_sums='passed'))


def validate_first_sample_against_full_original_dump(sid,panels):
    checkpoint=OUT/'first_sample_original_dump_validation.json'
    if checkpoint.exists():
        assert json.loads(checkpoint.read_text())['sample_id']==sid
        return
    observed={}
    with gzip.open(OUT/(sid+'.regional_raw_counts.tsv.gz'),'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):observed[r['kmer']]=int(r['retained_read_count'])
    expected={};source_rows=0
    for g in 'FM':
        dest=WORK/sid/('original_'+g+'.verification_dump.tsv')
        tool([TOOLS,'-t2','-hp','transform',INPUT/sid/(g+'_candidate_counts'),'dump',dest],WORK/sid/(g+'.independent_dump.log'))
        with dest.open() as f:
            for line in f:
                k,n=line.split();source_rows+=1
                if k in panels:
                    assert panels[k]['candidate_group']==g
                    expected[k]=int(n)
        remove_scratch(dest)
    assert expected==observed
    check_input(sid)
    save(checkpoint,dict(sample_id=sid,completed_at=now(),original_candidate_rows_streamed=source_rows,regional_nonzero_counts_independently_matched=len(observed),validation='Exact raw counts from full original dumps equal panel intersection output; source checksums unchanged'))


def main():
    WORK.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    lock=WORK/'runner.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    (WORK.parent/'regional_sensitivity_runner.pid').write_text(str(os.getpid()))
    try:
        verify_binary(KMC, 'kmc_sha256')
        plan=dict(regions=REGIONS,thresholds=THRESHOLDS,individuals=184,sex_counts={'F':92,'M':92},panel_selection='All unique exact-reference matches among fixed discovery-consistent Bonferroni heldout associations; start inside predefined region',
                  data_source='READ ONLY retained per-individual F/M candidate databases; no original reads counted again',
                  count_mode='KMC intersect preserves RIGHT input individual raw count; query panel input ci1; source input ci2',
                  validation='Input SHA256 before and after; panel canonical sequences/count1 checked; regional threshold2 prevalence sums against prior mapping',
                  inference='Descriptive threshold sensitivity and within-sex depth/QC correlations; no new sex-association P values, no local-depth-adjusted claim')
        plan=json.loads(json.dumps(plan))
        if (OUT/'analysis_plan.json').exists():assert json.loads((OUT/'analysis_plan.json').read_text())==plan
        else:save(OUT/'analysis_plan.json',plan)
        state('preparing_fixed_regional_panels')
        panels,metadata=prepare_panels()
        with (ROOT/'outputs/kmer/heldout_validation_manifest_184.tsv').open(encoding='utf-8-sig') as f:samples=list(csv.DictReader(f,delimiter='\t'))
        assert len(samples)==184 and collections.Counter(s['sex'] for s in samples)=={'F':92,'M':92}
        results=[]
        for s in samples:
            state('reading_retained_candidate_counts',current_sample=s['sample_id'],completed=len(results),total=184)
            results.append(sample(s,panels,metadata))
            if len(results)==1:validate_first_sample_against_full_original_dump(s['sample_id'],panels)
            print(now(),s['sample_id'],'regional counts complete',len(results),'/184',flush=True)
        state('independent_prevalence_validation_and_summary',completed=184,total=184)
        summarize(results,metadata)
        state('regional_threshold_sensitivity_complete_pending_biological_interpretation',completed=184,total=184)
    except Exception as exc:
        state('failed_preserved_for_checkpoint_inspection',error=repr(exc));raise
    finally:
        if lock.read_text().strip()==str(os.getpid()):lock.unlink()


if __name__=='__main__':main()
