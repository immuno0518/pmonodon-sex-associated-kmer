from release_config import *
"""Prespecified descriptive screen in 8F/8M; held-out validation is required.

Counts in prevalence databases represent individuals with read count >=2,
not sequencing depth. Retain in >=6/8 focal-sex samples and <=1/8 opposite-sex
samples. This is an association candidate screen, never a causal/W assignment.
"""
import csv,json,pathlib,subprocess,time,shutil
ROOT = RUN_ROOT
OUT=ROOT/'outputs/kmer'
RUN=ROOT/'work/kmer'
TOOL=KMC_TOOLS_BIN

def run(args,log):
    with log.open('w') as f:subprocess.run([str(TOOL),'-t2','-hp']+[str(x) for x in args],stdout=f,stderr=subprocess.STDOUT,check=True)

def delete_db(prefix):
    prefix=pathlib.Path(prefix).resolve()
    assert prefix.parent==(RUN/'discovery').resolve()
    for suffix in ['.kmc_pre','.kmc_suf']:
        p=pathlib.Path(str(prefix)+suffix)
        if p.exists():p.unlink()

def main():
    d=RUN/'discovery';d.mkdir(exist_ok=True)
    rows=list(csv.DictReader((OUT/'discovery_manifest_16.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
    assert all((RUN/'qc_full'/r['sample_id']/'complete.json').exists() for r in rows)
    if shutil.disk_usage(d).free<120*2**30:raise RuntimeError('Need 120 GiB free for prevalence operations')
    dbs={}
    for sex in ['F','M']:
        subset=[r for r in rows if r['sex']==sex];assert len(subset)==8
        prev=None
        for i,r in enumerate(subset):
            binary=d/(sex+'_binary_'+str(i));target=d/(sex+'_prevalence_'+str(i))
            source=RUN/'qc_full'/r['sample_id']/'k31'
            run(['transform',source,'-ci2','set_counts','1',binary],d/(sex+'_binary_'+str(i)+'.log'))
            if prev is None:prev=binary;continue
            run(['simple',prev,binary,'union',target,'-ocsum','-cs255'],d/(sex+'_union_'+str(i)+'.log'))
            delete_db(prev);delete_db(binary);prev=target
        dbs[sex]=prev
        run(['transform',prev,'histogram',OUT/(sex+'_discovery_prevalence_histogram.tsv')],d/(sex+'_histogram.log'))
    for sex,opposite in [('F','M'),('M','F')]:
        candidates=d/(sex+'_associated_candidates')
        run(['simple',dbs[sex],'-ci6',dbs[opposite],'-ci2','kmers_subtract',candidates],d/(sex+'_candidate_screen.log'))
        run(['transform',candidates,'dump',OUT/(sex+'_candidate_kmers_discovery.tsv')],d/(sex+'_dump.log'))
    (OUT/'discovery_screen_status.json').write_text(json.dumps(dict(phase='complete',F=8,M=8,k=31,
        per_sample_min_read_count=2,min_focal_prevalence=6,max_opposite_prevalence=1,
        count_column='number of focal-sex individuals with count >=2',
        heldout_validation='184 other individuals reserved; not yet tested',
        limitations='Discovery screen only. No statistical significance, genomic mapping, W assignment, or universal sex specificity inferred.'),indent=2),encoding='utf-8')
if __name__=='__main__':main()
