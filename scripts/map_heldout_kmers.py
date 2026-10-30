from release_config import *
"""Local exact, both-strand mapping of all replicated k-mers to the full reference.

Bounded query chunks keep memory moderate. Each completed chunk is checkpointed;
no original reads, candidate databases, or reference files are modified/uploaded.
The first chunk prioritizes perfect held-out separation; every other replicated
k-mer is mapped afterward, without duplication or changing association tests.
"""
import argparse
import collections
import csv
import gc
import gzip
import hashlib
import json
import math
import os
import pathlib
import sys
import time

ROOT = RUN_ROOT
import ahocorasick
import psutil

STATS = ROOT / 'outputs/kmer/heldout_statistics'
RESULT = ROOT / 'outputs/kmer/heldout_mapping'
WORK = ROOT / 'work/kmer/heldout_mapping'
STATUS = ROOT / 'outputs/kmer/heldout_mapping_status.json'
REF = input_path('reference')
REFERENCE_MANIFEST = ROOT / 'inputs/reference_sequence_manifest.tsv'
CHUNK_SIZE = 250000
QUERY_COLUMNS = ['query_id','candidate_group','kmer','heldout_F_present','heldout_M_present','bonferroni_p','perfect_heldout_separation']
MAP_COLUMNS = QUERY_COLUMNS + ['exact_reference_hits','unique_chr','unique_start_1based','unique_end_1based','unique_strand','first_3_hits','mononucleotide_entropy_bits']
RC = str.maketrans('ACGT', 'TGCA')


def digest(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def atomic_json(path,value):
    path=pathlib.Path(path); temp=path.with_suffix(path.suffix+'.partial')
    temp.write_text(json.dumps(value,indent=2),encoding='utf-8'); temp.replace(path)


def status(phase,**values):
    atomic_json(STATUS,dict(phase=phase,updated_at=time.strftime('%Y-%m-%d %H:%M:%S'),**values))


def read_fasta(path):
    name=None;parts=[]
    with pathlib.Path(path).open('r',encoding='ascii') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name,''.join(parts)
                name=line[1:].split()[0];parts=[]
            else: parts.append(line.strip())
    if name is not None: yield name,''.join(parts)


def prepare_queries():
    source=STATS/'direction_consistent_bonferroni_associations.tsv.gz'
    expected=json.loads((STATS/'result_checksums.json').read_text())[source.name]
    assert digest(source)==expected
    manifest_path=WORK/'query_manifest.json'
    if manifest_path.exists():
        manifest=json.loads(manifest_path.read_text())
        assert manifest['source_sha256']==expected
        for item in manifest['chunks']:
            assert digest(WORK/item['file'])==item['sha256']
        return manifest
    perfect_path=WORK/'queries_000_perfect.tsv'
    perfect_tmp=pathlib.Path(str(perfect_path)+'.partial')
    regular_files=[];regular_handle=None;regular_writer=None;regular_n=0
    counts=collections.Counter(); perfect_counts=collections.Counter(); n=0
    with perfect_tmp.open('w',newline='',encoding='ascii') as pf, gzip.open(source,'rt') as f:
        pw=csv.DictWriter(pf,fieldnames=QUERY_COLUMNS,delimiter='\t');pw.writeheader()
        for row in csv.DictReader(f,delimiter='\t'):
            group=row['candidate_group']; counts[group]+=1;n+=1
            nf,nm=int(row['heldout_nF']),int(row['heldout_nM'])
            fp,mp=int(row['heldout_F_present']),int(row['heldout_M_present'])
            perfect=(fp==nf and mp==0) if group=='F' else (mp==nm and fp==0)
            assert float(row['bonferroni_p_all_candidates'])<=.05 and row['direction_consistent']=='1'
            q=dict(query_id=group+'_%07d'%counts[group],candidate_group=group,kmer=row['kmer'],
                   heldout_F_present=fp,heldout_M_present=mp,
                   bonferroni_p=row['bonferroni_p_all_candidates'],perfect_heldout_separation=int(perfect))
            if perfect:
                pw.writerow(q);perfect_counts[group]+=1
            else:
                if regular_n % CHUNK_SIZE==0:
                    if regular_handle: regular_handle.close()
                    path=WORK/('queries_%03d.tsv'%(len(regular_files)+1))
                    regular_files.append(path)
                    regular_handle=pathlib.Path(str(path)+'.partial').open('w',newline='',encoding='ascii')
                    regular_writer=csv.DictWriter(regular_handle,fieldnames=QUERY_COLUMNS,delimiter='\t');regular_writer.writeheader()
                regular_writer.writerow(q);regular_n+=1
            if n%500000==0: status('preparing_queries',processed=n)
    if regular_handle:regular_handle.close()
    perfect_tmp.replace(perfect_path)
    for path in regular_files:pathlib.Path(str(path)+'.partial').replace(path)
    summary=json.loads((STATS/'association_summary.json').read_text())
    for group in ('F','M'):
        assert counts[group]==summary['candidate_groups'][group]['direction_consistent_bonferroni_significant']
        assert perfect_counts[group]==summary['candidate_groups'][group]['perfect_separation_in_heldout']
    items=[]
    for path in [perfect_path]+regular_files:
        with path.open() as f: number=sum(1 for _ in f)-1
        items.append(dict(file=path.name,count=number,sha256=digest(path)))
    manifest=dict(source_sha256=expected,total_queries=n,counts=dict(counts),perfect_counts=dict(perfect_counts),chunks=items)
    atomic_json(manifest_path,manifest)
    return manifest


def reference_identity():
    rows=list(csv.DictReader(REFERENCE_MANIFEST.open(encoding='utf-8-sig'),delimiter='\t'))
    expected={r['original_id']:r for r in rows}
    assert len(expected)==expected_sequences() and sum(int(r['length_bp']) for r in rows)==expected_reference_bp()
    return expected


def make_automaton(queries):
    automaton=ahocorasick.Automaton()
    for i,q in enumerate(queries):
        key=q['kmer'];assert len(key)==31 and set(key)<=set('ACGT')
        reverse=key.translate(RC)[::-1]
        assert key!=reverse and key not in automaton and reverse not in automaton
        automaton.add_word(key,2*i);automaton.add_word(reverse,2*i+1)
    automaton.make_automaton()
    return automaton


def map_chunk(item,index,total,expected):
    checkpoint=WORK/('chunk_%03d.complete.json'%index)
    output=RESULT/('chunk_%03d.mapping.tsv.gz'%index)
    identity=dict(query_sha256=item['sha256'],reference_map_sha256=digest(REFERENCE_MANIFEST),algorithm='exact_31mer_both_strands_v1')
    if checkpoint.exists():
        saved=json.loads(checkpoint.read_text());assert saved['identity']==identity
        assert digest(output)==saved['mapping_sha256']
        return saved
    queries=list(csv.DictReader((WORK/item['file']).open(),delimiter='\t'))
    assert len(queries)==item['count']
    status('building_exact_match_index',chunk=index+1,total_chunks=total,queries=len(queries))
    automaton=make_automaton(queries)
    assert psutil.Process().memory_info().rss < 4*2**30
    hits=[0]*len(queries);examples=[[] for _ in queries]
    by_chrom={};seen=set();bp=0
    print('Mapping chunk',index+1,'/',total,'queries',len(queries),'index',automaton.get_stats(),flush=True)
    for chrom,seq in read_fasta(REF):
        assert chrom in expected and chrom not in seen;seen.add(chrom)
        assert len(seq)==int(expected[chrom]['length_bp'])
        assert hashlib.sha256(seq.encode('ascii')).hexdigest()==expected[chrom]['sequence_sha256']
        text=seq.upper();groups=collections.Counter()
        for end,code in automaton.iter(text):
            i,reverse=divmod(code,2)
            hits[i]+=1;groups[queries[i]['candidate_group']]+=1
            if len(examples[i])<3:examples[i].append((chrom,end-29,end+1,'-' if reverse else '+'))
        bp+=len(seq);by_chrom[chrom]=dict(groups)
        status('mapping_full_reference',chunk=index+1,total_chunks=total,chromosome=chrom,
               scanned_bp=bp,total_reference_bp=expected_reference_bp(),queries=len(queries))
        print('Chunk',index+1,chrom,'verified and scanned',len(seq),'bp',flush=True)
    assert seen==set(expected) and bp==expected_reference_bp()
    categories={group:collections.Counter() for group in ('F','M')}
    partial=pathlib.Path(str(output)+'.partial')
    with gzip.open(partial,'wt',compresslevel=1,newline='') as out:
        writer=csv.DictWriter(out,fieldnames=MAP_COLUMNS,delimiter='\t');writer.writeheader()
        for i,q in enumerate(queries):
            category='unmapped' if hits[i]==0 else ('unique' if hits[i]==1 else 'multiple')
            categories[q['candidate_group']][category]+=1
            entropy=-sum((count/31)*math.log2(count/31) for count in collections.Counter(q['kmer']).values())
            unique=examples[i][0] if hits[i]==1 else ('','','','')
            writer.writerow(dict(q,exact_reference_hits=hits[i],unique_chr=unique[0],unique_start_1based=unique[1],
                                 unique_end_1based=unique[2],unique_strand=unique[3],
                                 first_3_hits=';'.join('%s:%d-%d:%s'%x for x in examples[i]),
                                 mononucleotide_entropy_bits='%.6f'%entropy))
    partial.replace(output)
    saved=dict(identity=identity,query_count=len(queries),mapping_sha256=digest(output),
               categories={g:dict(v) for g,v in categories.items()},reference_hits_by_chromosome=by_chrom,
               reference_sequences_verified=len(seen),reference_bp=bp)
    atomic_json(checkpoint,saved)
    del automaton,queries,hits,examples;gc.collect()
    return saved


def merge_results(manifest):
    status('summarizing_genomic_localization')
    reference_lengths={name:int(row['length_bp']) for name,row in reference_identity().items()}
    categories=collections.defaultdict(collections.Counter);chrom_counts=collections.Counter()
    perfect_categories=collections.defaultdict(collections.Counter);sdr=collections.Counter()
    intervals=collections.defaultdict(list);windows=collections.Counter();n=0
    target=RESULT/'all_replicated_candidates_reference_mapping.tsv.gz'
    partial=pathlib.Path(str(target)+'.partial')
    with gzip.open(partial,'wt',compresslevel=1,newline='') as out:
        writer=csv.DictWriter(out,fieldnames=MAP_COLUMNS,delimiter='\t');writer.writeheader()
        for i,item in enumerate(manifest['chunks']):
            with gzip.open(RESULT/('chunk_%03d.mapping.tsv.gz'%i),'rt') as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    writer.writerow(row);n+=1
                    group=row['candidate_group'];hits=int(row['exact_reference_hits'])
                    category='unmapped' if hits==0 else ('unique' if hits==1 else 'multiple')
                    categories[group][category]+=1
                    if row['perfect_heldout_separation']=='1':perfect_categories[group][category]+=1
                    if hits==1:
                        chrom=row['unique_chr'];start,end=int(row['unique_start_1based']),int(row['unique_end_1based'])
                        chrom_counts[(group,chrom)]+=1;intervals[(group,chrom)].append((start,end))
                        windows[(group,chrom,(start-1)//100000)]+=1
                        if chrom=='Chr19' and start<=2890000 and end>=2710000:sdr[group]+=1
    assert n==manifest['total_queries'];partial.replace(target)
    with (RESULT/'unique_mapping_chromosome_counts.tsv').open('w') as out:
        out.write('candidate_group\tchromosome\tunique_kmers\n')
        for (group,chrom),count in sorted(chrom_counts.items()):out.write('%s\t%s\t%d\n'%(group,chrom,count))
    with (RESULT/'unique_mapping_100kb_windows.tsv').open('w') as out:
        out.write('candidate_group\tchromosome\tstart_1based\tend_1based\tunique_kmers_by_start\n')
        for (group,chrom,window),count in sorted(windows.items()):
            out.write('%s\t%s\t%d\t%d\t%d\n'%(group,chrom,window*100000+1,min((window+1)*100000,reference_lengths[chrom]),count))
    with (RESULT/'unique_mapping_clusters_gap1kb.tsv').open('w') as out:
        out.write('candidate_group\tchromosome\tstart_1based\tend_1based\tspan_bp\tcovered_bp\tkmer_count\n')
        for (group,chrom),values in sorted(intervals.items()):
            values.sort();start=end=covered=number=0
            for left,right in values:
                if number and left>end+1000:
                    out.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\n'%(group,chrom,start,end,end-start+1,covered,number))
                    start=end=covered=number=0
                if not number:start=left;end=right;covered=right-left+1;number=1
                else:covered+=max(0,right-max(end,left-1));end=max(end,right);number+=1
            if number:out.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\n'%(group,chrom,start,end,end-start+1,covered,number))
    summary=dict(phase='exact_mapping_complete_pending_biological_interpretation',total_queries=n,
                 categories={g:dict(v) for g,v in categories.items()},
                 perfect_heldout_categories={g:dict(v) for g,v in perfect_categories.items()},
                 unique_kmers_overlapping_Chr19_2710000_2890000=dict(sdr),
                 reference_sequences=expected_sequences(),reference_bp=expected_reference_bp(),
                 mapping='All exact matches to both strands of the full local reference; case-insensitive; no mismatch allowance.',
                 limitations=['Unmapped does not imply W-specific sequence: reference gaps, haplotype differences, repeats, contamination or sequencing errors remain possible.',
                              'Adjacent and overlapping k-mers are correlated; k-mer counts are not counts of independent loci.',
                              'Single-family held-out association does not establish population-wide specificity or causality.',
                              'Per-k-mer hit counts include all exact reference occurrences; at most three example coordinates are displayed.',
                              'Unique mapping clusters join intervals separated by at most 1 kb; cluster spans can include uncovered gaps.'])
    atomic_json(RESULT/'mapping_summary.json',summary)
    atomic_json(RESULT/'final_mapping_checksum.json',{target.name:digest(target)})
    return summary


def self_test():
    queries=[{'kmer':'ACGTTGCAACGTAGCTAGGCTAACGTTACGA'}]
    assert len(queries[0]['kmer'])==31
    a=make_automaton(queries);key=queries[0]['kmer'];reverse=key.translate(RC)[::-1]
    text='NN'+key+'N'+reverse+'N'+key
    observed=[(end-29,end+1,'-' if code%2 else '+') for end,code in a.iter(text)]
    assert observed==[(3,33,'+'),(35,65,'-'),(67,97,'+')],observed
    assert not list(a.iter('N'*100))
    repetitive=make_automaton([{'kmer':'A'*31}])
    assert [(end-29,code%2) for end,code in repetitive.iter('A'*32+'N'+'T'*31)]==[(1,0),(2,0),(34,1)]
    print('Exact matching, both strands, repeated matches and 1-based coordinates verified.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--self-test',action='store_true');args=ap.parse_args()
    if args.self_test:self_test();return
    WORK.mkdir(exist_ok=True);RESULT.mkdir(exist_ok=True)
    lock=WORK/'runner.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    try:
        expected=reference_identity();manifest=prepare_queries()
        for i,item in enumerate(manifest['chunks']):map_chunk(item,i,len(manifest['chunks']),expected)
        summary=merge_results(manifest)
        status(summary['phase'],summary=str(RESULT/'mapping_summary.json'))
        print(json.dumps(summary),flush=True)
    except Exception as exc:
        status('failed',error=str(exc));raise
    finally:lock.unlink()


if __name__=='__main__':main()
