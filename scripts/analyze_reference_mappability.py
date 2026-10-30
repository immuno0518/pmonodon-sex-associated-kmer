from release_config import *
"""Local-only validation and descriptive opportunity-normalized regional analysis."""
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
import numpy as np

BASE = ROOT / 'outputs/kmer/reference_mappability'
DOWNLOAD = BASE / 'server_download'
MAP = ROOT / 'outputs/kmer/heldout_mapping'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n', encoding='utf8')


def main():
    plan = {
        'created_at': datetime.now().astimezone().isoformat(),
        'scope': 'Local-only descriptive rates of associated unique k31 starts among all reference-unique k31 starts, k31/e0/both strands',
        'regions': ['Chr19 whole', 'Chr19:2710000-2890000', 'Chr34 whole', 'Contig6 whole', 'all reference sequences'],
        'coordinates': '1-based inclusive; count kmer start in region to match opportunity denominator',
        'windows': 'Fixed nonoverlapping 100kb bins from coordinate 1; zero candidate bins included',
        'descriptive_window_rank_eligibility': 'Full 100000bp windows with at least 10000 reference-unique k31 starts',
        'no_regional_enrichment_pvalues': 'Overlapping k-mers, linked loci and posthoc region selection preclude an independent-binomial interpretation',
        'checks': ['All downloaded file SHA256 and sizes', 'Decompressed frequency SHA256', 'All 93 sequence lengths and sequence hashes against input reference manifest', 'Every 100kb opportunity count and every chromosome total recomputed from raw frequency', 'Every unique candidate maps to frequency1 and every reported multiple-hit example agrees with saturated frequency'],
        'unresolved_confounding': ['Single full-sib family', 'Sequencing batch not recorded', 'Overall read yield does not establish local depth equivalence'],
    }
    plan_path = BASE / 'local_analysis_plan.json'
    if plan_path.exists():
        old = json.loads(plan_path.read_text(encoding='utf8'))
        assert old['descriptive_window_rank_eligibility'] == plan['descriptive_window_rank_eligibility']
    else:
        save(plan_path, plan)
    manifest = json.loads((DOWNLOAD/'download_manifest.json').read_text())
    for r in manifest['files']:
        p = DOWNLOAD / Path(r['path']).name
        assert p.stat().st_size == r['bytes'] and sha(p) == r['sha256'], str(p)
    local_names = identity_names()
    seqs = json.loads((DOWNLOAD/'reference_sequence_manifest.json').read_text())['sequences']
    assert len(seqs) == len(local_names) == expected_sequences()
    starts, lengths, neutral_to_local = {}, {}, {}
    offset = 0
    for r in seqs:
        loc = local_names[r['sequence']]
        assert r['length_bp'] == int(loc['length_bp']) and r['sequence_sha256'] == loc['sequence_sha256']
        c = loc['original_id']
        starts[c], lengths[c] = offset, r['length_bp']
        neutral_to_local[r['sequence']] = c
        offset += r['length_bp']
    assert offset == expected_reference_bp()
    raw = BASE/'reference_k31_e0.freq8'
    if raw.exists():
        assert raw.stat().st_size == offset and sha(raw) == manifest['raw_frequency_sha256']
    else:
        temp = raw.with_suffix('.freq8.part')
        h = hashlib.sha256()
        with gzip.open(DOWNLOAD/'data01_k31_e0.freq8.gz', 'rb') as source, temp.open('wb') as dest:
            for b in iter(lambda: source.read(8*1024*1024), b''):
                h.update(b)
                dest.write(b)
        assert temp.stat().st_size == offset and h.hexdigest() == manifest['raw_frequency_sha256']
        temp.replace(raw)
    freq = np.memmap(str(raw), dtype=np.uint8, mode='r')
    opportunities, windows = {}, {}
    with (DOWNLOAD/'reference_k31_opportunity_100kb.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            c = neutral_to_local[r['sequence']]
            a, b = int(r['start_1based']), int(r['end_1based'])
            data = freq[starts[c]+a-1:starts[c]+b]
            valid, unique = int(np.count_nonzero(data)), int(np.count_nonzero(data == 1))
            assert valid == int(r['valid_k31_starts']) and unique == int(r['unique_k31_starts'])
            assert valid-unique == int(r['multiple_k31_starts'])
            assert int(np.count_nonzero(data == 255)) == int(r['saturated_frequency_positions'])
            windows[(c, (a-1)//100000)] = {'chromosome':c, 'start_1based':a, 'end_1based':b, 'length_bp':b-a+1, 'valid_k31_starts':valid, 'unique_k31_starts':unique, 'F_unique':0, 'M_unique':0, 'F_perfect_heldout_unique':0, 'M_perfect_heldout_unique':0}
    for r in json.loads((DOWNLOAD/'reference_k31_opportunity_summary.json').read_text())['sequences']:
        c = neutral_to_local[r['sequence']]
        data = freq[starts[c]:starts[c]+lengths[c]]
        valid, unique = int(np.count_nonzero(data)), int(np.count_nonzero(data==1))
        assert r['valid_k31_starts'] == valid and r['unique_k31_starts'] == unique
        opportunities[c] = {'valid_k31_starts':valid, 'unique_k31_starts':unique}
    counts = defaultdict(Counter)
    target = Counter()
    unique_checked, multi_checked, target_intervals = 0, 0, defaultdict(list)
    with gzip.open(MAP/'all_replicated_candidates_reference_mapping.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            hits = int(r['exact_reference_hits'])
            g = r['candidate_group']
            perfect = r['perfect_heldout_separation'] == '1'
            if hits == 1:
                c, a = r['unique_chr'], int(r['unique_start_1based'])
                assert freq[starts[c]+a-1] == 1
                unique_checked += 1
                counts[c][g] += 1
                counts[c][g+'_perfect'] += int(perfect)
                w = windows[(c, (a-1)//100000)]
                w[g+'_unique'] += 1
                w[g+'_perfect_heldout_unique'] += int(perfect)
                if c == 'Chr19' and 2710000 <= a <= 2890000:
                    target[g] += 1
                    target[g+'_perfect'] += int(perfect)
                    target_intervals[g].append((a, min(a+30, 2890000)))
            elif hits > 1:
                for hit in r['first_3_hits'].split(';'):
                    c, span, strand = hit.split(':')
                    a = int(span.split('-')[0])
                    assert freq[starts[c]+a-1] == min(255,hits), (hit,hits)
                    multi_checked += 1
    summary = json.loads((MAP/'mapping_summary.json').read_text())
    for g in 'FM':
        assert sum(r[g] for r in counts.values()) == summary['categories'][g]['unique']
    genome_unique = sum(r['unique_k31_starts'] for r in opportunities.values())
    genome_counts = {g:sum(r[g] for r in counts.values()) for g in 'FM'}
    baseline = {g:genome_counts[g]/genome_unique for g in 'FM'}

    def row(c, a, b, opp, cnt):
        r = dict(chromosome=c, start_1based=a, end_1based=b, length_bp=b-a+1, **opp)
        for g in 'FM':
            r[g+'_unique'] = cnt[g]
            r[g+'_perfect_heldout_unique'] = cnt[g+'_perfect']
            r[g+'_associated_per_million_unique_starts'] = cnt[g]/opp['unique_k31_starts']*1e6 if opp['unique_k31_starts'] else None
            r[g+'_density_fold_vs_genome'] = cnt[g]/opp['unique_k31_starts']/baseline[g] if opp['unique_k31_starts'] else None
        return r

    chrom_rows = [row(c,1,lengths[c],opportunities[c],counts[c]) for c in lengths]
    data = freq[starts['Chr19']+2710000-1:starts['Chr19']+2890000]
    target_row = row('Chr19',2710000,2890000,{'valid_k31_starts':int(np.count_nonzero(data)), 'unique_k31_starts':int(np.count_nonzero(data==1))},target)
    for g in 'FM':
        prev, union = 0, 0
        for a,b in sorted(target_intervals[g]):
            union += max(0,b-max(a-1,prev))
            prev = max(prev,b)
        target_row[g+'_union_covered_bp'] = union
    eligible = []
    for w in windows.values():
        for g in 'FM':
            w[g+'_associated_per_million_unique_starts'] = w[g+'_unique']/w['unique_k31_starts']*1e6 if w['unique_k31_starts'] else None
        if w['length_bp']==100000 and w['unique_k31_starts']>=10000:
            eligible.append(w)
    with (BASE/'unique_opportunity_normalized_100kb_windows.tsv').open('w',newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(next(iter(windows.values())).keys()), delimiter='\t')
        writer.writeheader()
        writer.writerows(windows.values())
    with (BASE/'unique_opportunity_normalized_chromosomes.tsv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(chrom_rows[0].keys()),delimiter='\t')
        writer.writeheader()
        writer.writerows(chrom_rows)
    result = dict(created_at=datetime.now().astimezone().isoformat(), reference_bp=offset,
                  reference_unique_k31_starts=genome_unique, total_associated_unique=genome_counts,
                  genome_associated_per_million_unique_starts={g:baseline[g]*1e6 for g in 'FM'},
                  target=target_row, chromosomes=chrom_rows,
                  eligible_100kb_windows=len(eligible),
                  top_20_normalized_100kb_windows={g:sorted(eligible,key=lambda w:w[g+'_associated_per_million_unique_starts'],reverse=True)[:20] for g in 'FM'},
                  interpretation='Opportunity-normalized spatial density is descriptive. No independent-locus regional enrichment P value is asserted.',
                  limitations=plan['unresolved_confounding']+['Exact reference uniqueness does not establish read alignment accuracy in alternate haplotypes.', 'Unmapped k-mers do not establish W linkage.', 'Perfect sample-split separation is within one family, not independent-population validation.'])
    save(BASE/'opportunity_normalized_regional_assessment.json',result)
    audit=dict(checked_at=datetime.now().astimezone().isoformat(), downloaded_checksums='passed', decompressed_sha256='passed', local_reference_sequence_hashes_and_lengths=str(expected_sequences()) + '/' + str(expected_sequences()) + ' passed', all_100kb_and_chromosome_opportunity_totals='passed', unique_candidate_frequency_checks=unique_checked, multiple_hit_example_frequency_checks=multi_checked, candidate_counts_match_full_mapping='passed', identity_sequence_names_used=True, server_data_deleted=False)
    save(BASE/'independent_local_validation.json',audit)
    print(json.dumps({'validation':audit,'reference_unique_k31_starts':genome_unique,'target':target_row,'requested_chromosomes':[r for r in chrom_rows if r['chromosome'] in ('Chr19','Chr34','Contig6')]},indent=2))


if __name__=='__main__':
    main()
