from release_config import *
from reviewer_genomics_common import *
import numpy as np, math
LD=OUT/'local_LD';GD=OUT/'gene_evidence'
data={g:np.load(LD/(g+'_DP5_pairwise_r2.npz')) for g in ('all200','female100','male100')}
common=sorted(set.intersection(*(set(d['positions']) for d in data.values())))
mat={}
for g,d in data.items():
    ix=np.searchsorted(d['positions'],common);mat[g]=d['r2'][np.ix_(ix,ix)]
valid=np.triu(np.ones((len(common),len(common)),bool),1)
for m in mat.values():valid &=np.isfinite(m)
matched={'common_informative_sites':len(common),'matched_valid_pairs':int(valid.sum()),'groups':{g:{'mean_r2':float(m[valid].mean()),'median_r2':float(np.median(m[valid])),'fraction_ge_0_8':float((m[valid]>=.8).mean())} for g,m in mat.items()},'interpretation':'Same SNP pairs in all groups, descriptive only. Different-sample r2 sampling behavior and within-family structure remain.'}
save(LD/'matched_pair_group_comparison.json',matched)
# Exact conditional 2x3 genotype test. Enumerate all tables with fixed margins.
def exact23(a,b):
    n=sum(a);cols=[x+y for x,y in zip(a,b)];N=sum(cols)
    observed=math.prod(math.comb(c,x) for c,x in zip(cols,a));den=math.comb(N,n)
    total=0;tail=0
    for x in range(max(0,n-cols[1]-cols[2]),min(cols[0],n)+1):
        for y in range(max(0,n-x-cols[2]),min(cols[1],n-x)+1):
            z=n-x-y;w=math.comb(cols[0],x)*math.comb(cols[1],y)*math.comb(cols[2],z);total+=w
            if w<=observed:tail+=w
    assert total==den
    return tail/den
rows=list(csv.DictReader((GD/'coding_splice_variants_all200.tsv').open(),delimiter='\t'))
n=sum(r['association_QC_eligible']=='True' for r in rows)
for r in rows:
    if r['association_QC_eligible']=='True':
        a=[int(r['F_'+k]) for k in ('ref_hom','het','alt_hom')];b=[int(r['M_'+k]) for k in ('ref_hom','het','alt_hom')]
        p=exact23(a,b);r['genotype_2x3_exact_P']=p;r['genotype_2x3_Bonferroni_P']=min(1,p*n)
    else:r['genotype_2x3_exact_P']='NA';r['genotype_2x3_Bonferroni_P']='NA'
tsv(GD/'coding_splice_variants_all200.tsv',rows)
tsv(GD/'predicted_high_impact_candidates.tsv',[r for r in rows if r['effect'].startswith('predicted_') or r['effect']=='CDS_boundary_complex'])
save(GD/'genotype_association_validation.json',{'method':'Two-sided exact conditional 2x3 genotype table test; all tables with fixed margins; tables with probability <= observed included; exact integer Vandermonde sum checked for every tested table','tests':n,'Bonferroni_scope':'all QC-eligible coding/splice records in these two preselected genes; not genome-wide','warning':'Single full-sib family, sex-dependent mapping/callability; association does not assign W/Z or establish function.'})
print(json.dumps(matched,indent=2))
