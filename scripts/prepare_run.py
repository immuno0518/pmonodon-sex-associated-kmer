"""Prepare directories and validate user-supplied manifests; copy no sequence data."""
import argparse, csv, hashlib, json, shutil
from release_config import *

def read_samples(path):
    with path.open(encoding='utf-8-sig') as f:r=list(csv.DictReader(f,delimiter='\t'))
    if not r:raise ValueError('Empty sample manifest')
    if not {'sample_id','sex','R1','R2','sequencing_batch'} <= set(r[0]):raise ValueError('Missing sample columns')
    ids=[x['sample_id'] for x in r]
    if len(set(ids))!=len(ids) or any(not x.replace('-','').isalnum() for x in ids):raise ValueError('Invalid or duplicate sample ID')
    if any(x['sex'] not in ('F','M') for x in r):raise ValueError('sex must be explicitly provided as F or M')
    for x in r:
        for k in ('R1','R2'):
            x[k]=str(resolve_input(x[k]))
            if not Path(x[k]).is_file():raise FileNotFoundError(x[k])
    return r

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifests',action='store_true');ap.add_argument('--reference',action='store_true');a=ap.parse_args()
    for name in ['inputs','work/kmer/discovery','work/kmer/qc_full','work/kmer/qc_pilot','work/kmer/heldout_validation','work/kmer/heldout_statistics','work/kmer/heldout_mapping','work/kmer/regional_sensitivity','work/matplotlib_cache','outputs/kmer/heldout_statistics','outputs/kmer/heldout_mapping','outputs/kmer/reference_mappability','outputs/kmer/regional_sensitivity','outputs/kmer/manuscript_main_figure','outputs/kmer/manuscript_main_figure_v2','outputs/kmer/manuscript_main_figure_v3','outputs/reviewer_followup_20261021']:
        (RUN_ROOT/name).mkdir(parents=True,exist_ok=True)
    if a.manifests:
        from collections import Counter
        d=read_samples(input_path('discovery_manifest'));h=read_samples(input_path('heldout_manifest'));allrows=read_samples(input_path('all_samples_manifest'))
        assert len(d)==16 and Counter(x['sex'] for x in d)=={'F':8,'M':8}
        assert len(h)==184 and Counter(x['sex'] for x in h)=={'F':92,'M':92}
        assert not {x['sample_id'] for x in d}&{x['sample_id'] for x in h}
        assert len(allrows)==200 and Counter(x['sex'] for x in allrows)=={'F':100,'M':100}
        assert {x['sample_id'] for x in allrows}=={x['sample_id'] for x in d+h}
        byid={x['sample_id']:x for x in d+h}
        assert all(x==byid[x['sample_id']] for x in allrows),'Inconsistent labels or read paths'
        for name,rows in [('discovery_manifest_16.tsv',d),('heldout_validation_manifest_184.tsv',h),('sample_manifest_200.tsv',allrows)]:
            p=RUN_ROOT/'outputs/kmer'/name
            assert not p.exists(),'Preserve existing manifest: '+str(p)
            with p.open('w',newline='',encoding='utf-8') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    if a.reference:
        from map_heldout_kmers import read_fasta
        rows=[]
        for name,seq in read_fasta(input_path('reference')):
            rows.append(dict(original_id=name,length_bp=len(seq),sequence_sha256=hashlib.sha256(seq.encode('ascii')).hexdigest()))
        assert len({r['original_id'] for r in rows})==len(rows)==expected_sequences()
        assert sum(r['length_bp'] for r in rows)==expected_reference_bp()
        dest=RUN_ROOT/'inputs/reference_sequence_manifest.tsv';assert not dest.exists()
        with dest.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    print('Preparation complete; inputs remain read-only.')
if __name__=='__main__':main()
