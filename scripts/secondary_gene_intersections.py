from pathlib import Path
import csv,json,hashlib
from release_config import *
KM=RUN_ROOT/'outputs/kmer'
OUT=RUN_ROOT/'outputs/annotation'
OUT.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(n,x):(OUT/n).write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def rows(p):
    with Path(p).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
clusters = rows(KM/'heldout_mapping/unique_mapping_clusters_gap1kb.tsv')
selected = []
for chrom in ['Chr34','Contig6']:
    for direction in ['F','M']:
        candidates = [r for r in clusters if r['chromosome']==chrom and r['candidate_group']==direction]
        if candidates:
            selected.append(max(candidates,key=lambda r:int(r['kmer_count'])))
GFF = input_path('gff')
genes = []
with GFF.open(encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):continue
        c=line.rstrip('\n').split('\t')
        if len(c)<9 or c[2]!='gene' or c[0] not in ['Chr34','Contig6']:continue
        attr=dict(kv.split('=',1) for kv in c[8].split(';') if '=' in kv)
        for peak in selected:
            if peak['chromosome']==c[0] and int(c[3])<=int(peak['end_1based']) and int(c[4])>=int(peak['start_1based']):
                genes.append({'gene_id':attr.get('ID',''),'chromosome':c[0], 'gene_start':int(c[3]),'gene_end':int(c[4]),
                              'strand':c[6], 'candidate_group':peak['candidate_group'],
                              'cluster_start':int(peak['start_1based']),'cluster_end':int(peak['end_1based']),
                              'cluster_kmers':int(peak['kmer_count'])})
ANNO = input_path('functional_annotations')
want={g['gene_id'] for g in genes}|set(CONFIG['annotation_focal_ids'])
ann={}
with ANNO.open(encoding='utf-8-sig', newline='') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['Gene_id'] in want and r['Gene_id'] not in ann:ann[r['Gene_id']]=r
for r in genes:
    a=ann.get(r['gene_id'],{})
    for k in ['NR_Accession','NR_Description','SwissProt_Gene_Name','Pfam_IDs','Pfam_Descriptions']:r[k]=a.get(k,'not available')
dump('secondary_cluster_gene_annotation.json',{'selection':'Highest candidate-count cluster per sex-association direction and secondary sequence; descriptive gene-body intersection only, not an independent candidate-gene test.',
     'clusters':selected,'overlapping_genes':genes,'GFF_sha256':sha(GFF),'annotation_sha256':sha(ANNO),
     'limits':'A missing gene-body intersection does not exclude nearby regulatory elements. Model IDs and inferred functional labels do not establish sex-determining function.'})
focal=[]
for gid in sorted(want):
    if gid in CONFIG['annotation_focal_ids']:
        r=ann.get(gid,{})
        focal.append({k:r.get(k,'not available') for k in ['Gene_id','Chrom','Gene_start','Gene_end','Strand','Query_length','NR_Accession','NR_Description','SwissProt_Gene_Name','Pfam_IDs','Pfam_Descriptions']})
dump('focal_gene_annotations.json',focal)
