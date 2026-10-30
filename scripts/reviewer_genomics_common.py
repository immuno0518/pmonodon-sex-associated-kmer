from release_config import *
"""Local read-only indexed VCF access; no individual data leave this machine."""
import sys, pathlib, gzip, struct, hashlib, json, csv
ROOT = RUN_ROOT
from Bio import bgzf
OUT=ROOT/'outputs/reviewer_followup_20261021'
SOURCE=input_path('variant_source_directory')
OLD=input_path('previous_CDS_directory')
REF=input_path('reference')
GFF=input_path('gff')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def save(p,x):
    p=pathlib.Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def tsv(p,rows,fields=None):
    p=pathlib.Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def bins(beg,end):
    end-=1
    return [0]+[i for off,shift in [(1,26),(9,23),(73,20),(585,17),(4681,14)] for i in range(off+(beg>>shift),off+(end>>shift)+1)]
def chunks(path,chrom,beg,end):
    with gzip.open(str(path)+'.tbi','rb') as f:
        assert f.read(4)==b'TBI\x01'
        n,fmt,cs,cb,ce,meta,skip,ln=struct.unpack('<8i',f.read(32))
        names=f.read(ln).rstrip(b'\0').decode().split('\0'); wanted=set(bins(beg,end)); result=[]
        for name in names:
            nb,=struct.unpack('<i',f.read(4))
            for _ in range(nb):
                bn,nc=struct.unpack('<Ii',f.read(8)); raw=f.read(nc*16)
                if name==chrom and bn in wanted: result += list(struct.iter_unpack('<QQ',raw))
            ni,=struct.unpack('<i',f.read(4));f.read(ni*8)
    merged=[]
    for a,b in sorted(result):
        if merged and a<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],b)
        else: merged.append([a,b])
    return merged
def indexed_lines(path,chrom,beg,end):
    # beg/end are zero-based half-open; return records whose POS is inside.
    seen=set()
    with bgzf.BgzfReader(str(path),'rb') as f:
        for a,b in chunks(path,chrom,beg,end):
            f.seek(a)
            while f.tell()<b:
                off=f.tell();line=f.readline()
                if not line:break
                if line.startswith(b'#'):continue
                v=line.split(b'\t',2)
                if v[0].decode()!=chrom:continue
                pos=int(v[1])
                if beg<pos<=end and off not in seen:
                    seen.add(off);yield line.decode()
def header(path):
    out=[]
    with bgzf.BgzfReader(str(path),'rb') as f:
        while True:
            line=f.readline().decode()
            if not line.startswith('#'):break
            out.append(line)
            if line.startswith('#CHROM'):break
    return out
def reference(chrom,beg,end):
    entries={x.split('\t')[0]:x.split('\t') for x in pathlib.Path(str(REF)+'.fai').read_text().splitlines()}
    x=entries[chrom]; offset,lb,lw=map(int,x[2:5]);buf=[]
    with REF.open('rb') as f:
        p=beg
        while p<end:
            n=min(end-p,lb-p%lb);f.seek(offset+(p//lb)*lw+p%lb);buf.append(f.read(n));p+=n
    return b''.join(buf).decode().upper()
def genotypes(v,min_dp=5,min_gq=20):
    fmt=v[8].split(':');idx={k:i for i,k in enumerate(fmt)};out=[]
    for raw in v[9:]:
        a=raw.strip().split(':');gt=a[idx['GT']].replace('|','/').split('/')
        ok=len(gt)==2 and all(x in ('0','1') for x in gt)
        for field,limit in [('DP',min_dp),('GQ',min_gq)]:
            if field in idx:
                z=a[idx[field]] if idx[field]<len(a) else '.'
                ok=ok and z!='.' and float(z)>=limit
        out.append(sum(map(int,gt)) if ok else -1)
    return out
if __name__=='__main__':
    OUT.mkdir(exist_ok=True,parents=True)
    plan={'region':'Chr19:2500001-3100000','SDR':'Chr19:2710000-2890000 (1-based inclusive)','samples':'all 200, F100/M100; no rediscovery or kmer rerun',
    'LD':'Unphased diploid allele-count correlation squared, pairwise-complete observations, pooled and separate sexes; genotype DP>=5 and GQ>=20 where fields supplied; sensitivity DP>=3; biallelic SNPs only; per-sex call rate>=0.90, pooled MAF>=0.05; within-sex analyses require within-sex MAF>=0.05 and nonzero genotype variance; pair completeness>=90%; no HWE filter or LD pruning; all eligible pairs aggregated in fixed 10kb coordinate bins and 10kb distance bins',
    'LD_limits':'Single full-sib family; descriptive sample genotype correlations, not independent-population LD or recombination boundaries; pooled sex contrasts may inflate r2; monomorphic/constant heterozygote genotypes undefined, not zero.',
    'gene26':'Audit annotated CDS/splice variants from existing local VCFs and available phased evidence; no Z/W label for unphased exons; predicted consequence is not functional proof.',
    'duplication_selection':'Independently reconstruct reference CDS, translate, align homologous protein segments and back-translate codons. Report descriptive homologous-region divergence only when valid; branch-specific selection requires reliable ortholog/paralog sequences and phylogeny. No inferred TE-mediated mechanism before verified annotation/breakpoints.',
    'privacy':'Individual data are user-supplied inputs; all source files read-only; new outputs in this dedicated directory.'}
    save(OUT/'analysis_plan.json',plan)
    manifest=[]
    for kind in ('snp','indel'):
        src=SOURCE/('06.pop_filter_vcf/pop_filter.'+kind+'.vcf.gz');dest=OUT/('Chr19_2500001_3100000.'+kind+'.vcf.gz')
        if dest.exists():raise RuntimeError('Refusing to overwrite extracted VCF '+str(dest))
        hh=header(src);count=0;formats=set()
        with gzip.open(dest,'wt',encoding='utf-8') as f:
            f.writelines(hh)
            for line in indexed_lines(src,'Chr19',2500000,3100000):
                f.write(line);count+=1;formats.add(line.split('\t')[8])
        # Different, wider chunk request must recover exactly the same focal records.
        direct=[x for x in indexed_lines(src,'Chr19',2490000,3110000) if 2500000<int(x.split('\t')[1])<=3100000]
        with gzip.open(dest,'rt') as f: saved=[x for x in f if not x.startswith('#')]
        assert direct==saved
        manifest.append({'kind':kind,'source':str(src),'source_size':src.stat().st_size,'source_mtime_ns':src.stat().st_mtime_ns,'index_sha256':sha(str(src)+'.tbi'),'output':str(dest),'output_sha256':sha(dest),'records':count,'sample_count':len(hh[-1].rstrip().split('\t'))-9,'formats':sorted(formats),'wider_index_query_exact_match':True})
        print(kind,count,formats,flush=True)
    save(OUT/'indexed_extraction_validation.json',manifest)
