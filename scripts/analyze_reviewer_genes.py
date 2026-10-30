from release_config import *
from reviewer_genomics_common import *
import numpy as np, collections
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio.codonalign.codonseq import CodonSeq,cal_dn_ds
from scipy.stats import fisher_exact
GOUT=OUT/'gene_evidence';GOUT.mkdir(exist_ok=True,parents=True)
DOUT=OUT/'duplication_selection';DOUT.mkdir(exist_ok=True,parents=True)
genes={}
for line in GFF.open():
    if line.startswith('#'):continue
    v=line.rstrip().split('\t')
    if len(v)<9 or v[0]!='Chr19':continue
    for label in ('gene24','gene26'):
        gid=CONFIG['focal_gene_ids'][label]
        if gid not in v[8]:continue
        if v[2]=='gene':genes[label]={'id':gid,'start':int(v[3]),'end':int(v[4]),'strand':v[6],'cds_exons':[]}
        if v[2]=='CDS':genes[label]['cds_exons'].append((int(v[3]),int(v[4]),int(v[7])))
refcheck=[]
for label,g in genes.items():
    ex=sorted(g['cds_exons'],reverse=g['strand']=='-');g['cds_exons']=ex
    coords=[p for a,b,phase in ex for p in (range(a,b+1) if g['strand']=='+' else range(b,a-1,-1))]
    g['coords']=coords;g['map']={p:i for i,p in enumerate(coords)}
    g['cds']=''.join(reference('Chr19',a-1,b) if g['strand']=='+' else str(Seq(reference('Chr19',a-1,b)).reverse_complement()) for a,b,phase in ex)
    # Subsequent exon phase describes split codons; do not trim each exon.
    assert ex[0][2]==0 and len(g['cds'])%3==0
    previous=str(SeqIO.read(OLD/('02_haplotype_sequences/'+label+'_reference.CDS.fa'),'fasta').seq)
    assert g['cds']==previous,(label,'reference reconstruction mismatch')
    protein=str(Seq(g['cds']).translate());g['protein']=protein.rstrip('*')
    assert '*' not in g['protein']
    (DOUT/(label+'_reference_reconstructed.CDS.fa')).write_text('>'+label+'\n'+g['cds']+'\n')
    refcheck.append({'gene':label,'CDS_length':len(g['cds']),'protein_length_without_terminal_stop':len(g['protein']),'matches_existing_CDS':True,'terminal_stop':protein.endswith('*'),'internal_stop_count':g['protein'].count('*'),'strand':g['strand'],'exon_count':len(ex)})
samplesex={r['sample_id']:r['sex'] for r in csv.DictReader((ROOT/'outputs/kmer/sample_manifest_200.tsv').open(encoding='utf-8-sig'),delimiter='\t')}
def consequence(v,g):
    p=int(v[1]);ref,alt=v[3:5];end=p+len(ref)-1
    if ',' in alt or any(x not in 'ACGT' for x in ref+alt):return 'complex_or_multiallelic_not_annotated','NA','NA'
    actual=reference('Chr19',p-1,end)
    assert actual==ref,(p,ref,actual)
    idx=[g['map'].get(q) for q in range(p,end+1)]
    fully=all(i is not None for i in idx)
    if fully:
        lo=min(idx);hi=max(idx)+1
        rr=ref if g['strand']=='+' else str(Seq(ref).reverse_complement());aa=alt if g['strand']=='+' else str(Seq(alt).reverse_complement())
        assert g['cds'][lo:hi]==rr
        mut=g['cds'][:lo]+aa+g['cds'][hi:];delta=len(alt)-len(ref)
        if delta%3:return 'predicted_frameshift',lo+1,'NA'
        prot=str(Seq(mut[:len(mut)//3*3]).translate());base=str(Seq(g['cds']).translate())
        if base[0]=='M' and prot[0]!='M':effect='predicted_start_lost'
        elif '*' in prot[:-1] and ('*' not in base[:-1] or prot.index('*')<base.index('*')):effect='predicted_stop_gained'
        elif base.endswith('*') and not prot.endswith('*'):effect='predicted_stop_lost'
        elif delta:effect='inframe_indel'
        elif prot==base:effect='synonymous'
        else:effect='missense'
        aa_change=('%s%d%s'%(base[lo//3],lo//3+1,prot[lo//3])) if len(ref)==len(alt)==1 else 'NA'
        return effect,lo+1,aa_change
    # Canonical two intronic bases at each internal exon junction, strand aware.
    ex=sorted(g['cds_exons'])
    for (a,b,_),(c,d,__) in zip(ex,ex[1:]):
        for x,y,side in [(b+1,b+2,'donor' if g['strand']=='+' else 'acceptor'),(c-2,c-1,'acceptor' if g['strand']=='+' else 'donor')]:
            if p<=y and end>=x:return 'predicted_splice_'+side,'NA','NA'
    if any(i is not None for i in idx):return 'CDS_boundary_complex','NA','NA'
    return 'intronic','NA','NA'
records=[];allbody=collections.Counter();highimpact=[]
for kind in ('snp','indel'):
    with gzip.open(OUT/('Chr19_2500001_3100000.'+kind+'.vcf.gz'),'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:];continue
            if line.startswith('#'):continue
            v=line.rstrip().split('\t');p=int(v[1])
            for label,g in genes.items():
                if not(g['start']<=p<=g['end']):continue
                effect,cpos,aachange=consequence(v,g);allbody[label]+=1
                if effect=='intronic':continue
                row={'gene':label,'position':p,'ref':v[3],'alt':v[4],'FILTER':v[6],'effect':effect,'CDS_position':cpos,'AA_change_single_variant':aachange,'interpretation':'single-variant reference-model prediction; not a phased Z/W allele or functional assay'}
                if ',' not in v[4]:
                    gt=genotypes(v);n=np.array(gt)
                    for s in ('F','M'):
                        y=n[[i for i,k in enumerate(samples) if samplesex[k]==s]];valid=y>=0
                        row[s+'_called']=int(valid.sum());row[s+'_ref_hom']=int((y==0).sum());row[s+'_het']=int((y==1).sum());row[s+'_alt_hom']=int((y==2).sum());row[s+'_carrier']=int((y>0).sum());row[s+'_call_rate']=float(valid.mean())
                    a,b=row['F_carrier'],row['M_carrier'];nf,nm=row['F_called'],row['M_called']
                    row['carrier_Fisher_P']=float(fisher_exact([[a,nf-a],[b,nm-b]])[1]) if nf and nm else None
                    row['association_QC_eligible']=nf>=90 and nm>=90 and v[6] in ('.','PASS')
                else:
                    for s in ('F','M'):
                        for k in ('called','ref_hom','het','alt_hom','carrier','call_rate'):row[s+'_'+k]='NA'
                    row['carrier_Fisher_P']=None;row['association_QC_eligible']=False
                records.append(row)
                if effect.startswith('predicted_') or effect=='CDS_boundary_complex':highimpact.append(row)
n_tests=sum(r['association_QC_eligible'] for r in records)
for r in records:r['carrier_Bonferroni_P']=min(1,r['carrier_Fisher_P']*n_tests) if r['association_QC_eligible'] else 'NA'
tsv(GOUT/'coding_splice_variants_all200.tsv',records)
if highimpact:tsv(GOUT/'predicted_high_impact_candidates.tsv',highimpact)
phasefile=input_path('phased_variants')
phase=list(csv.DictReader(phasefile.open(),delimiter='\t'));phasecoding=[]
for r in phase:
    for label,g in genes.items():
        if int(r['pos']) in g['map']:phasecoding.append(dict(gene=label,**r))
if phasecoding:tsv(GOUT/'existing_phased_CDS_variants.tsv',phasecoding)
candidate_matches=[r for r in records if r['position']==2737771 and r['ref']=='T' and r['alt']=='TG']
save(GOUT/'summary.json',{'status':'population_coding_splice_screen_completed','reference_reconstruction':refcheck,'gene_body_variant_records':dict(allbody),'coding_splice_records':len(records),'effect_counts':{g:dict(collections.Counter(r['effect'] for r in records if r['gene']==g)) for g in genes},'association_tests_QC_eligible':n_tests,'phased_CDS_variant_counts':dict(collections.Counter(r['gene'] for r in phasecoding)),'gene24_2737771_T_TG_in_population_filtered_VCF':bool(candidate_matches),'Z_W_gene26_CDS_resolved':False,'limitations':['Existing SNP/indel callset ascertainment; absence of a call is not absence of the variant.','Single-variant consequences do not establish complete haplotype ORFs.','Missense does not imply deleterious; predicted LOF requires read/splice validation.','Within-family sex association is not independent W linkage or a functional demonstration.']})
aligner=Align.PairwiseAligner();aligner.mode='local';aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
p1,p2=genes['gene24']['protein'],genes['gene26']['protein'];al=aligner.align(p1,p2)[0];inds=al.indices
(DOUT/'reference_protein_local_alignment.txt').write_text(str(al),encoding='utf-8')
columns=[]
for k,(i,j) in enumerate(inds.T):
    columns.append({'alignment_column':k+1,'gene24_aa_position':int(i+1) if i>=0 else 'gap','gene26_aa_position':int(j+1) if j>=0 else 'gap','gene24_aa':p1[i] if i>=0 else '-','gene26_aa':p2[j] if j>=0 else '-','gene24_codon':genes['gene24']['cds'][3*i:3*i+3] if i>=0 else '---','gene26_codon':genes['gene26']['cds'][3*j:3*j+3] if j>=0 else '---'})
tsv(DOUT/'protein_guided_codon_alignment.tsv',columns)
estimates=[]
for trim in (0,2,5):
    gapcols=[k for k,(i,j) in enumerate(inds.T) if i<0 or j<0]
    kept=[k for k,(i,j) in enumerate(inds.T) if i>=0 and j>=0 and all(abs(k-h)>trim for h in gapcols)]
    c1=''.join(columns[k]['gene24_codon'] for k in kept);c2=''.join(columns[k]['gene26_codon'] for k in kept)
    assert all(str(Seq(columns[k]['gene24_codon']).translate())==columns[k]['gene24_aa'] and str(Seq(columns[k]['gene26_codon']).translate())==columns[k]['gene26_aa'] for k in kept)
    for name,seq in [('gene24',c1),('gene26',c2)]:
        with (DOUT/('homologous_codons_gapflank%d.fa'%trim)).open('w' if name=='gene24' else 'a') as f:f.write('>'+name+'\n'+seq+'\n')
    dn,ds=cal_dn_ds(CodonSeq(c1),CodonSeq(c2),method='NG86')
    estimates.append({'gap_flank_exclusion_aa':trim,'codons':len(kept),'dN_NG86':float(dn),'dS_NG86':float(ds),'omega':float(dn/ds) if ds>0 and dn>=0 else None,'interpretation':'Exploratory pairwise divergence of alignable reference paralog segments; no branch/site selection test, no W/Z inference.'})
tsv(DOUT/'exploratory_pairwise_dN_dS.tsv',estimates)
non_gap=(inds[0]>=0)&(inds[1]>=0);ident=sum(p1[i]==p2[j] for i,j in inds.T if i>=0 and j>=0)
save(DOUT/'summary.json',{'status':'reference_homology_and_exploratory_pairwise_divergence_completed','reference_CDS_validation':refcheck,'genomic_gap_between_genes_bp':genes['gene26']['start']-genes['gene24']['end']-1,'gene24_strand':genes['gene24']['strand'],'gene26_strand':genes['gene26']['strand'],'protein_alignment_columns':len(columns),'ungapped_homologous_residue_pairs':int(non_gap.sum()),'identical_amino_acids':ident,'identity_among_residue_pairs':ident/int(non_gap.sum()),'identity_including_gap_columns':ident/len(columns),'gene24_aa_span':[int(inds[0][inds[0]>=0].min()+1),int(inds[0].max()+1)],'gene26_aa_span':[int(inds[1][inds[1]>=0].min()+1),int(inds[1].max()+1)],'exploratory_pairwise_estimates':estimates,'formal_branch_specific_selection_test':'not identifiable from only this reference paralog pair; reliable ortholog CDS and gene tree required','TE_duplication_mechanism':'pending verified full-reference TE annotation and breakpoint/homology evidence; proximity not causality','interpretation':'Local partial homology; genes on opposite strands with intervening gene25. Do not label a proven simple same-orientation tandem duplication or new regulatory function.'})
print(json.dumps({'coding_splice_records':len(records),'predicted_high_impact':len(highimpact),'dN_dS':estimates},indent=2),flush=True)
