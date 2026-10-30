from release_config import *
"""Independent read-back of durable regional raw counts and exportable summary."""
import collections
import csv
import gzip
import hashlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib_cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=ROOT/'outputs/kmer/regional_sensitivity'
WORK=ROOT/'work/kmer/regional_sensitivity'


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    status=json.loads((ROOT/'outputs/kmer/regional_sensitivity_status.json').read_text())
    assert status['phase']=='regional_threshold_sensitivity_complete_pending_biological_interpretation'
    with (OUT/'fixed_regional_panels.tsv').open() as f:panels={r['kmer']:r for r in csv.DictReader(f,delimiter='\t')}
    with (ROOT/'outputs/kmer/heldout_validation_manifest_184.tsv').open(encoding='utf-8-sig') as f:samples=list(csv.DictReader(f,delimiter='\t'))
    metadata=json.loads((WORK/'panels.complete.json').read_text())
    prevalence=collections.Counter();rows_checked=0;hashes={}
    for s in samples:
        sid=s['sample_id'];checkpoint=json.loads((WORK/sid/'complete.json').read_text())
        p=OUT/(sid+'.regional_raw_counts.tsv.gz')
        h=digest(p);assert h==checkpoint['raw_panel_counts_sha256'];hashes[p.name]=h
        hist=collections.defaultdict(collections.Counter);seen=set()
        with gzip.open(p,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                k=r['kmer'];assert k not in seen;seen.add(k)
                expected=panels[k]
                assert expected['region']==r['region'] and expected['candidate_group']==r['candidate_group']
                n=int(r['retained_read_count']);assert 2<=n<=65535
                hist[r['region']+'_'+r['candidate_group']][n]+=1
                rows_checked+=1
        for name,values in checkpoint['count_histograms'].items():
            assert dict(hist[name])=={int(k):v for k,v in values.items()}
        for row in checkpoint['summary_rows']:
            name=row['region']+'_'+row['candidate_group']
            observed=sum(v for n,v in hist[name].items() if n>=row['read_count_threshold'])
            assert observed==row['detected_kmers']
            if row['read_count_threshold']==2:prevalence[(name,s['sex'])]+=observed
    for name,x in metadata['expected'].items():
        for sex in 'FM':assert prevalence[(name,sex)]==x[sex+'_prevalence_sum']
    validation=dict(completed_at=datetime.now().astimezone().isoformat(),individual_compressed_files=184,raw_nonzero_count_rows=rows_checked,
                    all_compressed_sha256='passed',panel_membership_and_no_duplicates='passed',stored_count_histograms_and_all_threshold_rows='passed',
                    independent_mapped_prevalence_sums='passed',source_checksums_before_and_after='736 files passed in runner',raw_file_sha256=hashes)
    (OUT/'independent_compressed_counts_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    summary=json.loads((OUT/'regional_threshold_summary.json').read_text())
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axs=plt.subplots(2,3,figsize=(12,7),sharex=True,sharey=True)
    colors={'F':'#B63B83','M':'#2374AB'}
    names={'Chr19_target':'Chr19:2.71-2.89 Mb','Chr34':'Chr34','Contig6':'Contig6'}
    for col,region in enumerate(names):
        for row,g in enumerate('FM'):
            ax=axs[row,col];values=[r for r in summary['summary'] if r['region']==region and r['candidate_group']==g]
            for sex in 'FM':
                y=np.array([r[sex]['median'] for r in values])*100
                lo=np.array([r[sex]['q1'] for r in values])*100
                hi=np.array([r[sex]['q3'] for r in values])*100
                ax.plot(range(4),y,marker='o',markersize=4,color=colors[sex],label='Female individuals' if sex=='F' else 'Male individuals')
                ax.fill_between(range(4),lo,hi,color=colors[sex],alpha=.17)
            ax.set_title(names[region]+'\n'+format(values[0]['panel_kmers'],',')+' '+g+'-associated k-mers',fontweight='bold')
            ax.set_xticks(range(4));ax.set_xticklabels(['2','3','5','10'])
            ax.set_ylim(-2,102);ax.grid(axis='y',alpha=.2)
            if col==0:ax.set_ylabel('Panel detected per individual (%)')
            if row==1:ax.set_xlabel('Minimum retained read count')
    axs[0,0].legend(frameon=False,loc='center left',fontsize=8)
    fig.suptitle('Sensitivity of fixed regional k31 panels in 184 heldout individuals',fontsize=14,fontweight='bold',y=.995)
    fig.text(.015,.02,'Lines: median; bands: interquartile range (92 female, 92 male individuals). Raw count databases are reused; no reads recounted.\nPanels were selected after association and unique mapping. Threshold changes do not equalize depth or provide independent validation.\nSingle full-sib family; unrecorded batches; counts below 2 censored. No new association P-value or W-linkage claim.',fontsize=9)
    fig.tight_layout(rect=(0,.11,1,.965))
    for suffix in ('png','pdf'):fig.savefig(OUT/('regional_threshold_sensitivity.'+suffix),dpi=220)
    plt.close(fig)
    files=['analysis_plan.json','fixed_regional_panels.tsv','first_sample_original_dump_validation.json','per_individual_region_threshold_sensitivity.tsv','regional_threshold_summary.json','independent_prevalence_and_input_validation.json','independent_compressed_counts_validation.json','regional_threshold_sensitivity.png','regional_threshold_sensitivity.pdf']
    hashes.update({n:digest(OUT/n) for n in files})
    (OUT/'result_checksums.json').write_text(json.dumps({'created_at':datetime.now().astimezone().isoformat(),'sha256':hashes},indent=2)+'\n')
    compact=[{'region':r['region'],'candidate_group':r['candidate_group'],'threshold':r['read_count_threshold'],'F_median_fraction':r['F']['median'],'M_median_fraction':r['M']['median']} for r in summary['summary']]
    print(json.dumps({'validation_rows':rows_checked,'summary':compact},indent=2))


if __name__=='__main__':main()
