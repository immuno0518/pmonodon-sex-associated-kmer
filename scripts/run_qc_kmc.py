from release_config import *
"""Run validated native QC/KMC pipeline sequentially with space checks and logs.

Default: four complete paired-end samples, to measure full-sample resources.
This stage does not establish sex specificity. It does not launch the 200-sample cohort.
"""
import argparse,csv,json,os,pathlib,shutil,subprocess,sys,time,psutil
ROOT = RUN_ROOT
OUT=ROOT/'outputs/kmer'
RUN=ROOT/'work/kmer'
KMC=KMC_BIN
TOOLS=KMC_TOOLS_BIN
LIB=None

def logrun(cmd,path,env=None):
    peak_rss=0
    with path.open('w',encoding='utf-8') as log:
        p=subprocess.Popen([str(x) for x in cmd],stdout=log,stderr=subprocess.STDOUT,env=env)
        monitor=psutil.Process(p.pid)
        while p.poll() is None:
            try:
                rss=sum(x.memory_info().rss for x in [monitor]+monitor.children(recursive=True) if x.is_running())
                peak_rss=max(peak_rss,rss)
            except psutil.Error:pass
            time.sleep(0.5)
        code=p.returncode
    if code:raise RuntimeError('Command failed (%s): %s; see %s'%(code,cmd,path))
    return peak_rss

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['pilot','full'],default='pilot');ap.add_argument('--manifest',type=pathlib.Path,default=OUT/'pilot_manifest.tsv');args=ap.parse_args()
    rows=list(csv.DictReader(args.manifest.open(encoding='utf-8-sig'),delimiter='\t'))
    if args.stage=='full':
        unfinished=sum(not (RUN/'qc_full'/r['sample_id']/'complete.json').exists() for r in rows)
        # A measured first-sample database is approximately 10.2 GiB.
        if shutil.disk_usage(RUN).free < (120+unfinished*12)*2**30:
            raise RuntimeError('Requested retained databases exceed the measured disk budget')
    env=os.environ.copy();env['PYTHONIOENCODING']='utf-8'
    statuspath=OUT/'qc_kmc_run_status.json'
    results=[]
    for r in rows:
        sid=r['sample_id'];d=RUN/('qc_'+args.stage)/sid;d.mkdir(parents=True,exist_ok=True)
        tmp=d/'tmp';tmp.mkdir(exist_ok=True)
        # Keep ample room for a whole uncompressed paired sample and KMC spill files.
        reserve=120 if args.stage=='full' else 40
        if shutil.disk_usage(d).free<reserve*2**30:raise RuntimeError('Insufficient free disk space; keep %s GiB reserve'%reserve)
        inputs=[pathlib.Path(r[x]) if args.stage=='full' else RUN/'pilot'/sid/(x+'.fq') for x in ('R1','R2')]
        clean=[d/'R1.clean.fq',d/'R2.clean.fq']
        db=d/'k31';done=d/'complete.json'
        if done.exists():
            results.append(json.loads(done.read_text()));continue
        statuspath.write_text(json.dumps(dict(stage=args.stage,current_sample=sid,phase='QC',completed=[x['sample_id'] for x in results],started_at=time.strftime('%Y-%m-%d %H:%M:%S')),indent=2),encoding='utf-8')
        # MGI adapter prefixes independently observed in 197/200 company fastp reports.
        # Other three reports contained no adapter sequence >=20bp; a separate QC audit is retained.
        qc=[sys.executable,'-m','cutadapt','-j','2','-a','AAGTCGGAGGCCAAGCGGTC','-A','AAGTCGGATCGTAGCCATGT',
            '-O','8','-e','0.1','-q','20','--max-n','0.05','-m','50','--pair-filter=any',
            '--json',d/'cutadapt.json','-o',clean[0],'-p',clean[1],*inputs]
        start=time.time();qcpeak=logrun(qc,d/'cutadapt.log',env);qcsec=time.time()-start
        lst=d/'inputs.lst';lst.write_text('\n'.join(str(p) for p in clean)+'\n',encoding='utf-8')
        cmd=[KMC,'-k31','-m2','-sm','-t4','-ci2','-cs65535','-hp','-j'+str(d/'kmc.json'),'@'+str(lst),db,tmp]
        statuspath.write_text(json.dumps(dict(stage=args.stage,current_sample=sid,phase='KMC',completed=[x['sample_id'] for x in results]),indent=2),encoding='utf-8')
        start=time.time();kmcpeak=logrun(cmd,d/'kmc.log');kmcsec=time.time()-start
        logrun([TOOLS,'-t2','transform',db,'histogram',d/'k31.histogram.tsv'],d/'histogram.log')
        result=dict(sample_id=sid,sex=r['sex'],stage=args.stage,QC_seconds=qcsec,KMC_seconds=kmcsec,QC_peak_RSS_bytes=qcpeak,KMC_peak_RSS_bytes=kmcpeak,
            cutadapt_command=[str(x) for x in qc],kmc_command=[str(x) for x in cmd],
            database_bytes=sum(p.stat().st_size for p in d.glob('k31.kmc_*')),
            kmc_summary=json.loads((d/'kmc.json').read_text()),QC_summary=json.loads((d/'cutadapt.json').read_text()))
        done.write_text(json.dumps(result,indent=2),encoding='utf-8');results.append(result)
        (OUT/('qc_'+args.stage+'_metrics.json')).write_text(json.dumps(results,indent=2),encoding='utf-8')
        # Delete only the exact two task-generated transient files, after successful count.
        for p in clean:
            assert p.parent==d and p.name in ('R1.clean.fq','R2.clean.fq')
            p.unlink()
        print(sid,args.stage,'complete; QC',round(qcsec,1),'s; KMC',round(kmcsec,1),'s',flush=True)
    statuspath.write_text(json.dumps(dict(stage=args.stage,phase='complete',completed=[x['sample_id'] for x in results],note='Resource measurement only. No cohort association yet.'),indent=2),encoding='utf-8')
if __name__=='__main__':main()
