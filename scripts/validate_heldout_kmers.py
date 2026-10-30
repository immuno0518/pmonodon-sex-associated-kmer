from release_config import *
"""Sequential held-out validation; retain candidate counts, recycle full databases.

The discovery candidate set and count>=2 threshold are frozen before validation.
Absence from a candidate intersection means count<2, not proof of biological absence.
"""
import csv, hashlib, json, os, pathlib, shutil, subprocess, sys, time
from run_qc_kmc import logrun, KMC, TOOLS, LIB

ROOT = RUN_ROOT
OUT = ROOT / 'outputs/kmer'
RUN = ROOT / 'work/kmer'
BASE = RUN / 'heldout_validation'
STATUS = OUT / 'heldout_validation_status.json'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def atomic_json(path, value):
    temp = path.with_suffix('.json.partial')
    temp.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temp.replace(path)

def main():
    BASE.mkdir(exist_ok=True)
    # Exclusive process lock: remove only after checking an old process has exited.
    lock = BASE / 'runner.lock'
    with lock.open('x') as f: f.write(str(os.getpid()))
    try:
        rows = list(csv.DictReader((OUT/'heldout_validation_manifest_184.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
        assert len(rows)==184 and len({r['sample_id'] for r in rows})==184
        frozen = {}
        for sex in ('F','M'):
            p = OUT / (sex+'_candidate_kmers_discovery.tsv')
            with p.open() as f: count=sum(1 for _ in f)
            frozen[sex]={'sha256':digest(p),'candidate_count':count}
        freeze_path=BASE/'frozen_candidates.json'
        if freeze_path.exists(): assert json.loads(freeze_path.read_text())==frozen
        else: atomic_json(freeze_path,frozen)
        env=os.environ.copy();env['PYTHONIOENCODING']='utf-8'
        completed=[]
        for row in rows:
            sid=row['sample_id']; assert sid.replace('-','').isalnum()
            d=BASE/sid;d.mkdir(exist_ok=True);done=d/'complete.json'
            if done.exists(): completed.append(sid);continue
            if shutil.disk_usage(BASE).free<120*2**30: raise RuntimeError('Need 120 GiB temporary free space')
            def status(phase): atomic_json(STATUS,dict(phase=phase,current_sample=sid,completed=completed,total=184,updated_at=time.strftime('%Y-%m-%d %H:%M:%S')))
            status('QC')
            clean=[d/'R1.clean.fq',d/'R2.clean.fq'];tmp=d/'tmp';tmp.mkdir(exist_ok=True)
            summary_path=d/'count_summary.json';db=d/'k31'
            if not summary_path.exists():
                qc=[sys.executable,'-m','cutadapt','-j','2','-a','AAGTCGGAGGCCAAGCGGTC','-A','AAGTCGGATCGTAGCCATGT','-O','8','-e','0.1','-q','20','--max-n','0.05','-m','50','--pair-filter=any','--json',d/'cutadapt.json','-o',clean[0],'-p',clean[1],row['R1'],row['R2']]
                qc_checkpoint=d/'qc_complete.json'
                if qc_checkpoint.exists():
                    saved=json.loads(qc_checkpoint.read_text(encoding='utf-8'))
                    assert saved['cutadapt_command']==list(map(str,qc))
                    assert saved['clean_file_sizes']=={p.name:p.stat().st_size for p in clean}
                    assert saved['cutadapt_report_sha256']==digest(d/'cutadapt.json')
                    qcpeak=saved['QC_peak_RSS_bytes'];qcsec=saved['QC_seconds']
                else:
                    start=time.time();qcpeak=logrun(qc,d/'cutadapt.log',env);qcsec=time.time()-start
                    atomic_json(qc_checkpoint,dict(cutadapt_command=list(map(str,qc)),clean_file_sizes={p.name:p.stat().st_size for p in clean},cutadapt_report_sha256=digest(d/'cutadapt.json'),QC_peak_RSS_bytes=qcpeak,QC_seconds=qcsec))
                lst=d/'inputs.lst';lst.write_text('\n'.join(map(str,clean))+'\n',encoding='utf-8')
                cmd=[KMC,'-k31','-m2','-sm','-t4','-ci2','-cs65535','-hp','-j'+str(d/'kmc.json'),'@'+str(lst),db,tmp]
                status('KMC');start=time.time();peak=logrun(cmd,d/'kmc.log');elapsed=time.time()-start
                atomic_json(summary_path,dict(sample_id=sid,sex=row['sex'],QC_seconds=qcsec,KMC_seconds=elapsed,QC_peak_RSS_bytes=qcpeak,KMC_peak_RSS_bytes=peak,cutadapt_command=list(map(str,qc)),kmc_command=list(map(str,cmd))))
            status('candidate_intersections')
            checksums={}
            for sex in ('F','M'):
                target=d/(sex+'_candidate_counts')
                source=RUN/'discovery'/(sex+'_associated_candidates')
                cmd=[TOOLS,'-t2','-hp','simple',source,db,'-ci2','intersect',target,'-ocright','-cs65535']
                logrun(cmd,d/(sex+'_intersection.log'))
                hist=d/(sex+'_candidate_count_histogram.tsv')
                logrun([TOOLS,'-t2','transform',target,'histogram',hist],d/(sex+'_histogram.log'))
                for suffix in ('.kmc_pre','.kmc_suf'):
                    p=pathlib.Path(str(target)+suffix);assert p.is_file();checksums[p.name]=digest(p)
            atomic_json(done,dict(sample_id=sid,sex=row['sex'],candidate_database_sha256=checksums,frozen_candidates=frozen,per_sample_presence_threshold=2,absence_definition='read count <2 after QC; not biological absence'))
            # Reclaim only this worker's exact transient files after durable candidate outputs.
            for p in clean+[pathlib.Path(str(db)+s) for s in ('.kmc_pre','.kmc_suf')]:
                assert p.parent==d
                if p.exists(): p.unlink()
            completed.append(sid)
            print(sid,'candidate validation complete',flush=True)
        atomic_json(STATUS,dict(phase='candidate_counts_complete_pending_statistical_analysis',completed=completed,total=184,updated_at=time.strftime('%Y-%m-%d %H:%M:%S')))
    except Exception as exc:
        atomic_json(STATUS,dict(phase='failed',error=str(exc),updated_at=time.strftime('%Y-%m-%d %H:%M:%S')))
        raise
    finally:
        lock.unlink()

if __name__=='__main__': main()
