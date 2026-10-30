from release_config import *
"""Read existing discovery DBs for an all-200 descriptive figure; never recount reads."""
import csv, gzip, json, os, subprocess, hashlib
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
BASE=ROOT/'outputs/kmer'
OUT=BASE/'manuscript_main_figure_v3'
WORK=ROOT/'work/kmer/discovery_figure_panels'
TOOL=KMC_TOOLS_BIN
PANEL=ROOT/'work/kmer/regional_sensitivity'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def rows(p):
    with p.open(encoding='utf-8-sig') as f: return list(csv.DictReader(f,delimiter='\t'))

def save(p,obj):
    p.write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')

def run(args,log):
    with log.open('w') as f:
        subprocess.run([str(TOOL),'-t2','-hp']+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,check=True)

def fingerprint(prefix):
    result={}
    for ext in ('.kmc_pre','.kmc_suf'):
        p=Path(str(prefix)+ext);s=p.stat()
        result[ext]={'size':s.st_size,'mtime_ns':s.st_mtime_ns}
    return result

def main():
    WORK.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
    lock=WORK/'runner.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    try:
        verify_binary(TOOL, 'kmc_tools_sha256')
        meta=json.loads((PANEL/'panels.complete.json').read_text())
        for g in 'FM':
            for ext,digest in meta['panel_database_sha256'][g].items():assert sha(Path(str(PANEL/('panel_'+g))+ext))==digest
        panel_file=BASE/'regional_sensitivity/fixed_regional_panels.tsv'
        assert sha(panel_file)==meta['panel_file_sha256']
        panels={r['kmer']:r for r in rows(panel_file) if r['region']=='Chr19_target'}
        assert Counter(r['candidate_group'] for r in panels.values())=={'F':32753,'M':4507}
        discovery=rows(BASE/'discovery_manifest_16.tsv')
        heldout=rows(BASE/'heldout_validation_manifest_184.tsv')
        assert Counter(r['sex'] for r in discovery)=={'F':8,'M':8}
        assert not {r['sample_id'] for r in discovery}&{r['sample_id'] for r in heldout}
        save(OUT/'all200_display_plan.json',{'scope':'Add 16 discovery individuals to existing 184 heldout panel summaries; descriptive all-200 display, no new tests',
              'region_1based_inclusive':['Chr19',2710000,2890000],'read_count_threshold':2,'fixed_panel_sizes':{'F':32753,'M':4507},
              'discovery_per_sex':8,'heldout_per_sex':92,'source':'Existing qc_full/<sample>/k31 databases; no FASTQ read or counting run',
              'preservation_check':'Source size and mtime before/after read-only KMC intersection; not a new full-file SHA256 audit',
              'validation':'Every focal-sex k-mer prevalence compared against original discovery TSV; opposite-sex prevalence <=1',
              'limits':['Selected-panel descriptive summary','Discovery individuals are not validation individuals','Single family','Local depth effects remain']})
        union=WORK/'combined_regional_panel'
        if not (WORK/'union.complete.json').exists():
            run(['simple',PANEL/'panel_F',PANEL/'panel_M','union',union,'-ocmax'],WORK/'union.log')
            save(WORK/'union.complete.json',{'panel_sha256':meta['panel_file_sha256']})
        assert json.loads((WORK/'union.complete.json').read_text())['panel_sha256']==meta['panel_file_sha256']
        summaries=[];prevalence={'F':Counter(),'M':Counter()};source_records={}
        for i,s in enumerate(discovery):
            sid=s['sample_id'];sex=s['sex']; source=ROOT/'work/kmer/qc_full'/sid/'k31'
            complete=json.loads((source.parent/'complete.json').read_text())
            assert complete['sample_id']==sid and complete['sex']==sex and '-ci2' in complete['kmc_command'] and '-k31' in complete['kmc_command']
            before=fingerprint(source)
            d=WORK/sid;d.mkdir(exist_ok=True)
            dest=OUT/(sid+'.discovery_Chr19_raw_counts.tsv.gz');checkpoint=d/'complete.json'
            if not checkpoint.exists():
                intersect=d/'panel_intersection';dump=d/'panel_intersection.tsv'
                run(['simple',union,'-ci1',source,'-ci2','intersect',intersect,'-ocright','-cs65535'],d/'intersection.log')
                run(['transform',intersect,'dump',dump],d/'dump.log')
                with dump.open() as f,gzip.open(str(dest)+'.partial','wt') as out:
                    out.write('kmer\tcandidate_group\tretained_read_count\n')
                    for line in f:
                        k,n=line.split()
                        if k in panels:
                            assert 2<=int(n)<=65535
                            out.write(k+'\t'+panels[k]['candidate_group']+'\t'+n+'\n')
                assert fingerprint(source)==before
                Path(str(dest)+'.partial').replace(dest)
                save(checkpoint,{'source_fingerprint':before,'raw_counts_sha256':sha(dest),'panel_sha256':meta['panel_file_sha256']})
            ck=json.loads(checkpoint.read_text())
            assert ck['source_fingerprint']==before and ck['raw_counts_sha256']==sha(dest) and ck['panel_sha256']==meta['panel_file_sha256']
            detected=Counter();seen=set()
            with gzip.open(dest,'rt') as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    k=r['kmer'];g=r['candidate_group'];assert k not in seen and panels[k]['candidate_group']==g
                    assert int(r['retained_read_count'])>=2
                    seen.add(k);detected[g]+=1;prevalence[sex][k]+=1
            for g,n in [('F',32753),('M',4507)]:
                summaries.append({'sample_id':sid,'sex':sex,'cohort':'Discovery','candidate_group':g,'panel_kmers':n,
                                  'detected_kmers':detected[g],'detected_percent':detected[g]/n*100})
            source_records[sid]=ck
            print(datetime.now().isoformat(),i+1,'/16 discovery panels extracted',flush=True)
            save(WORK/'status.json',{'phase':'extracting_existing_counts','completed':i+1,'total':16,'pid':os.getpid()})
        focal_checks=0
        for g in 'FM':
            observed=set()
            with (BASE/(g+'_candidate_kmers_discovery.tsv')).open() as f:
                for line in f:
                    k,n=line.split()
                    if k in panels:
                        assert panels[k]['candidate_group']==g and prevalence[g][k]==int(n)
                        observed.add(k);focal_checks+=1
            assert observed=={k for k,r in panels.items() if r['candidate_group']==g}
            opposite='M' if g=='F' else 'F'
            assert all(prevalence[opposite][k]<=1 for k in observed)
        for r in rows(BASE/'regional_sensitivity/per_individual_region_threshold_sensitivity.tsv'):
            if r['region']=='Chr19_target' and r['read_count_threshold']=='2':
                summaries.append({'sample_id':r['sample_id'],'sex':r['sex'],'cohort':'Heldout','candidate_group':r['candidate_group'],
                                  'panel_kmers':int(r['panel_kmers']),'detected_kmers':int(r['detected_kmers']),
                                  'detected_percent':float(r['detected_fraction'])*100})
        assert len(summaries)==400
        assert Counter((r['candidate_group'],r['sex']) for r in summaries)=={(g,s):100 for g in 'FM' for s in 'FM'}
        with (OUT/'panels_DE_all200_detection.tsv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter='\t');writer.writeheader();writer.writerows(summaries)
        save(OUT/'discovery_panel_validation.json',{'status':'passed','discovery_individuals':16,'heldout_individuals':184,
             'focal_per_kmer_prevalences_independently_matched':focal_checks,'opposite_prevalence_le1_all_passed':True,
             'source_preservation':'Size/mtime unchanged; inputs used read-only; no full SHA audit of large discovery DBs',
             'sources_and_output_hashes':source_records,'all200_table_sha256':sha(OUT/'panels_DE_all200_detection.tsv')})
        save(WORK/'status.json',{'phase':'complete','completed':16,'total':16,'pid':os.getpid()})
    finally:
        if lock.read_text()==str(os.getpid()):lock.unlink()

if __name__=='__main__':main()
