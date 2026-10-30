from release_config import *
"""Synthetic-only GenMap validation; contains no study data or sequence names."""
import collections
import json
import random
import subprocess
from pathlib import Path

root = RUN_ROOT
tool = GENMAP_BIN
out = root / 'outputs/kmer/reference_mappability/server_download/toy_v2'
out.mkdir(parents=True, exist_ok=True)
random.seed(20261021)
seq = ''.join(random.choice('ACGT') for _ in range(150))
rc = lambda s: s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
seqs = {'test1': seq + 'N' * 31 + seq[:65], 'test2': rc(seq[45:120]) + 'N' * 20 + 'A' * 50,
        'test3': 'NNN' + seq[11:52] + 'NN', 'short': 'ACGTACGT', 'saturated_repeat': 'A' * 400}
fa = out / 'toy.fa'
fa.write_text(''.join('>' + c + '\n' + s + '\n' for c, s in seqs.items()))
subprocess.run([str(tool), 'index', '-F', str(fa), '-I', str(out / 'index'), '--version-check', '0'], check=True)
subprocess.run([str(tool), 'map', '-I', str(out / 'index'), '-O', str(out / 'toy'), '-K', '31', '-E', '0', '-T', '2', '-fs', '-t', '-bg', '-r', '--version-check', '0'], check=True)
freq = collections.Counter()
for s in seqs.values():
    for i in range(len(s) - 30):
        k = s[i:i+31]
        if set(k) <= set('ACGT'):
            freq[min(k, rc(k))] += 1
expected = {c: [min(255, freq[min(s[i:i+31], rc(s[i:i+31]))]) if len(s[i:i+31]) == 31 and set(s[i:i+31]) <= set('ACGT') else 0 for i in range(len(s))] for c,s in seqs.items()}
expected_bytes = bytes(v for c in seqs for v in expected[c])
assert (out / 'toy.freq8').read_bytes() == expected_bytes, 'raw frequency mismatch'
lines = (out / 'toy.txt').read_text().splitlines()
for i, c in enumerate(seqs):
    assert lines[i * 2] == '>' + c
    assert list(map(int, lines[i * 2 + 1].split())) == expected[c]
bg = {c: [0] * len(s) for c, s in seqs.items()}
for line in (out / 'toy.bedgraph').read_text().splitlines():
    if line.startswith(('track', '#')):
        continue
    c, a, b, v = line.split()
    a, b, v = int(a), int(b), int(v)
    bg[c][a:b] = [v] * (b-a)
assert bg == expected, 'bedgraph coordinates mismatch'
(out / 'expected.json').write_text(json.dumps(expected, indent=2))
result = {'passed': True, 'positions_checked': len(expected_bytes), 'checks': ['Forward and reverse-complement cross-contig counts', 'Overlapping repeat counts', 'Ambiguous-base windows and short-contig ends zero', 'Frequency8 saturates at 255 and cannot wrap to 1', 'Raw concatenation follows FASTA order', 'Bedgraph 0-based half-open coordinates', 'Raw, text, bedgraph agree with independent canonical Python counts']}
(out / 'validation.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
