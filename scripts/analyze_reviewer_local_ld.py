from release_config import *
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
os.environ['MKL_NUM_THREADS']='2'
from reviewer_genomics_common import *
import numpy as np
import collections, math, time
LD=OUT/'local_LD';LD.mkdir(exist_ok=True,parents=True)
rows=list(csv.DictReader((ROOT/'outputs/kmer/sample_manifest_200.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
sex={r['sample_id']:r['sex'] for r in rows}
variants=[];samples=None;reject=collections.Counter()
with gzip.open(OUT/'Chr19_2500001_3100000.snp.vcf.gz','rt') as f:
    for line in f:
        if line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:];continue
        if line.startswith('#'):continue
        v=line.rstrip().split('\t')
        if len(v[3])!=1 or len(v[4])!=1 or v[3] not in 'ACGT' or v[4] not in 'ACGT':reject['not_biallelic_ACGT_SNP']+=1;continue
        if v[6] not in ('.','PASS'):reject['FILTER_not_PASS']+=1;continue
        variants.append(v)
assert len(samples)==200 and set(samples)==set(sex)
groups={'all200':np.arange(200),'female100':np.array([i for i,s in enumerate(samples) if sex[s]=='F']),'male100':np.array([i for i,s in enumerate(samples) if sex[s]=='M'])}
pos=np.array([int(v[1]) for v in variants]);un,poscounts=np.unique(pos,return_counts=True);dups=set(un[poscounts>1]);dedup=np.array([p not in dups for p in pos])
summaries={};validation=[]
for dp in (5,3):
    geno=np.array([genotypes(v,dp,20) for v in variants],dtype=np.int8)
    mask=geno>=0;freq=np.maximum(geno,0).sum(1)/(2*np.maximum(mask.sum(1),1));maf=np.minimum(freq,1-freq)
    eligible=dedup & (maf>=.05)
    for name in ('female100','male100'):eligible &= mask[:,groups[name]].mean(1)>=.90
    site_rows=[]
    for i,v in enumerate(variants):
        site_rows.append(dict(position=int(pos[i]),ref=v[3],alt=v[4],call_rate_F=float(mask[i,groups['female100']].mean()),call_rate_M=float(mask[i,groups['male100']].mean()),pooled_MAF=float(maf[i]),eligible=bool(eligible[i])))
    tsv(LD/('site_QC_DP%d.tsv'%dp),site_rows)
    summaries[str(dp)]={'eligible_sites':int(eligible.sum()),'groups':{}}
    for group,inds in groups.items():
        x=geno[:,inds];m=x>=0;af=np.maximum(x,0).sum(1)/(2*np.maximum(m.sum(1),1));gmaf=np.minimum(af,1-af)
        keep=eligible & (gmaf>=.05)
        # A group can have constant heterozygous genotypes despite MAF=0.5.
        variance=np.array([np.var(a[a>=0]) if np.any(a>=0) else 0 for a in x]);constant=int((keep & (variance==0)).sum());keep &= variance>0
        x=x[keep].astype(float);pp=pos[keep];m=(x>=0).astype(float);x[x<0]=0;x2=x*x
        nsites=len(pp);bs=256;nb=60
        bsum=np.zeros(nb*nb);bcount=np.zeros(nb*nb,dtype=np.int64);bhigh=np.zeros(nb*nb,dtype=np.int64)
        dsum=np.zeros(nb);dcount=np.zeros(nb,dtype=np.int64);dhigh=np.zeros(nb,dtype=np.int64)
        raw=np.full((nsites,nsites),np.nan,dtype=np.float32) if dp==5 else None
        rng=np.random.RandomState(20261021);audit_budget=500;valid_pairs=0
        for a in range(0,nsites,bs):
            z=slice(a,min(a+bs,nsites));n=m[z]@m.T;sx=x[z]@m.T;sy=m[z]@x.T;sxx=x2[z]@m.T;syy=m[z]@x2.T;sxy=x[z]@x.T
            with np.errstate(divide='ignore',invalid='ignore'):
                cov=sxy-sx*sy/n;vx=sxx-sx*sx/n;vy=syy-sy*sy/n;r2=cov*cov/(vx*vy)
            valid=(n>=math.ceil(.9*len(inds)))&(vx>1e-10)&(vy>1e-10)&np.isfinite(r2)
            r2=np.clip(r2,0,1);r2[~valid]=np.nan
            if raw is not None:raw[z]=r2
            ii,jj=np.where(valid & (np.arange(a,a+len(x[z]))[:,None]<np.arange(nsites)[None,:]))
            gi=ii+a;vals=r2[ii,jj];valid_pairs+=len(vals)
            bi=(pp[gi]-2500001)//10000;bj=(pp[jj]-2500001)//10000;key=bi*nb+bj
            bcount+=np.bincount(key,minlength=nb*nb);bsum+=np.bincount(key,weights=vals,minlength=nb*nb);bhigh+=np.bincount(key,weights=(vals>=.8).astype(int),minlength=nb*nb).astype(np.int64)
            dist=(pp[jj]-pp[gi])//10000
            dcount+=np.bincount(dist,minlength=nb);dsum+=np.bincount(dist,weights=vals,minlength=nb);dhigh+=np.bincount(dist,weights=(vals>=.8).astype(int),minlength=nb).astype(np.int64)
            if dp==5 and len(vals) and audit_budget:
                take=rng.choice(len(vals),min(40,len(vals),audit_budget),replace=False)
                for k in take:
                    i,j=gi[k],jj[k];both=(m[i]>0)&(m[j]>0);independent=float(np.corrcoef(x[i,both],x[j,both])[0,1]**2)
                    err=abs(independent-float(vals[k]));assert err<1e-10,(group,i,j,err)
                    validation.append({'group':group,'pos1':int(pp[i]),'pos2':int(pp[j]),'paired_samples':int(both.sum()),'r2_moments':float(vals[k]),'r2_direct_corrcoef':independent,'absolute_error':err})
                audit_budget-=len(take)
        heat=[]
        for i in range(nb):
            for j in range(i,nb):
                k=i*nb+j;c=int(bcount[k]);heat.append({'bin1_start':2500001+i*10000,'bin2_start':2500001+j*10000,'pair_count':c,'mean_r2':bsum[k]/c if c else 'NA','fraction_r2_ge_0_8':bhigh[k]/c if c else 'NA'})
        decay=[{'distance_start_bp':i*10000,'distance_end_bp':(i+1)*10000-1,'pair_count':int(dcount[i]),'mean_r2':dsum[i]/dcount[i] if dcount[i] else 'NA','fraction_r2_ge_0_8':dhigh[i]/dcount[i] if dcount[i] else 'NA'} for i in range(nb)]
        tsv(LD/('%s_DP%d_heatmap_10kb.tsv'%(group,dp)),heat);tsv(LD/('%s_DP%d_distance_10kb.tsv'%(group,dp)),decay)
        if raw is not None:np.savez_compressed(LD/(group+'_DP5_pairwise_r2.npz'),positions=pp,r2=raw)
        summaries[str(dp)]['groups'][group]={'samples':len(inds),'informative_sites':nsites,'constant_genotype_sites_excluded':constant,'valid_pairs':valid_pairs,'mean_r2':float(dsum.sum()/dcount.sum()) if dcount.sum() else None}
        print('DP',dp,group,summaries[str(dp)]['groups'][group],flush=True)
tsv(LD/'independent_pairwise_validation.tsv',validation)
save(LD/'summary.json',{'status':'computed_validated','source_SNP_records':5782,'biallelic_PASS_records':len(variants),'rejected':dict(reject),'duplicate_position_sites_excluded':int((~dedup).sum()),'results':summaries,'independent_direct_pair_checks':len(validation),'max_absolute_error':max(x['absolute_error'] for x in validation),'note':'Descriptive unphased genotype r2 in one full-sib family. Undefined pairs remain NA; no recombination or population LD claim.'})
