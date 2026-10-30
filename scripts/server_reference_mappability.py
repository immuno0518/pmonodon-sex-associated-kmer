from release_config import *
"""Reference-only k31 opportunity analysis. No study individuals or local names."""
import csv
import gzip
import hashlib
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
OUT = ROOT / 'outputs/kmer/reference_mappability/server_download'
TOOL = GENMAP_BIN
REF = input_path('reference')
LOCK = OUT / 'runner.lock'
STATUS = OUT / 'status.json'


def stamp():
    return datetime.now().astimezone().isoformat()


def write_json(path, value):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def state(phase, **kwargs):
    write_json(STATUS, dict(checked_at=stamp(), phase=phase, pid=os.getpid(), **kwargs))
    print(stamp(), phase, flush=True)


def file_records(paths):
    return [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in paths]


def check_records(records):
    for r in records:
        p = ROOT / r['path']
        assert p.stat().st_size == r['bytes'] and sha(p) == r['sha256'], str(p)


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    # Atomic single-runner guard. Stale locks require explicit process verification.
    fd = os.open(str(LOCK), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    (OUT / 'runner.pid').write_text(str(os.getpid()))
    try:
        toy = json.loads((OUT / 'toy_v2/validation.json').read_text())
        assert toy['passed']
        verify_binary(TOOL, 'genmap_sha256')
        state('reference_manifest')
        records = []
        whole = hashlib.sha256()
        name = None
        length = 0
        digest = None
        with REF.open('rb') as f:
            for line in f:
                whole.update(line)
                if line.startswith(b'>'):
                    if name is not None:
                        records.append(dict(sequence=name, length_bp=length, sequence_sha256=digest.hexdigest()))
                    name = line[1:].split()[0].decode('ascii')
                    assert name, 'Reference sequence identifier must be nonempty'
                    length, digest = 0, hashlib.sha256()
                else:
                    seq = line.strip()
                    length += len(seq)
                    digest.update(seq)
        records.append(dict(sequence=name, length_bp=length, sequence_sha256=digest.hexdigest()))
        assert len(records) == expected_sequences() and sum(r['length_bp'] for r in records) == expected_reference_bp()
        reference_sha = whole.hexdigest()
        write_json(OUT / 'reference_sequence_manifest.json', dict(reference_sha256=reference_sha, sequences=records))
        plan = dict(created_at=stamp(), purpose='Exact k31 unique-start background from reference only; no individual or phenotype data',
                    genmap_version='1.3.0', genmap_binary_sha256=sha(TOOL), reference_sha256=reference_sha,
                    k=31, errors=0, both_strands=True, threads=8, niceness=10,
                    output='Frequency8 in original FASTA order; each sequence has length_bp bytes, 0 for invalid or terminal positions, 1 unique, 2..255 multiple (saturates at 255)',
                    window_coordinates='1-based inclusive windows; each value counts 31-mer starts in that window',
                    toy_validation=toy,
                    limitation='Reference unique-start opportunity is not a family/batch/depth-adjusted sex association model; overlapping k-mers are correlated.')
        plan_path = OUT / 'analysis_plan.json'
        if plan_path.exists():
            old = json.loads(plan_path.read_text())
            assert old['reference_sha256'] == reference_sha and old['genmap_binary_sha256'] == plan['genmap_binary_sha256']
        else:
            write_json(plan_path, plan)
        index = OUT / 'index'
        checkpoint = OUT / 'index.complete.json'
        if checkpoint.exists():
            check_records(json.loads(checkpoint.read_text())['files'])
        else:
            assert not index.exists(), 'Incomplete index preserved; inspect before resuming'
            state('indexing_reference', reference_bp=expected_reference_bp())
            subprocess.run([str(TOOL), 'index', '-F', str(REF), '-I', str(index), '-A', 'divsufsort', '--version-check', '0'], check=True)
            write_json(checkpoint, dict(completed_at=stamp(), reference_sha256=reference_sha, files=file_records(sorted(p for p in index.rglob('*') if p.is_file()))))
        raw = OUT / 'data01_k31_e0.freq8'
        checkpoint = OUT / 'mapping.complete.json'
        if checkpoint.exists():
            check_records(json.loads(checkpoint.read_text())['files'])
        else:
            assert not raw.exists(), 'Incomplete raw output preserved; inspect before resuming'
            state('computing_reference_k31_uniqueness', threads=8)
            subprocess.run([str(TOOL), 'map', '-I', str(index), '-O', str(OUT / 'data01_k31_e0'), '-K', '31', '-E', '0', '-T', '8', '-fs', '-r', '--version-check', '0'], check=True)
            assert raw.stat().st_size == sum(r['length_bp'] for r in records)
            write_json(checkpoint, dict(completed_at=stamp(), reference_sha256=reference_sha, files=file_records([raw])))
        state('summarizing_reference_opportunity')
        window_path = OUT / 'reference_k31_opportunity_100kb.tsv'
        seq_summary = []
        with raw.open('rb') as f, window_path.open('w') as out:
            out.write('sequence\tstart_1based\tend_1based\tvalid_k31_starts\tunique_k31_starts\tmultiple_k31_starts\tsaturated_frequency_positions\n')
            for r in records:
                n, c = r['length_bp'], r['sequence']
                s = dict(sequence=c, length_bp=n, valid_k31_starts=0, unique_k31_starts=0, multiple_k31_starts=0, saturated_frequency_positions=0)
                for a in range(0, n, 100000):
                    b = min(a+100000, n)
                    data = f.read(b-a)
                    assert len(data) == b-a
                    valid, unique, saturated = len(data)-data.count(b'\0'), data.count(b'\1'), data.count(b'\xff')
                    values = [valid, unique, valid-unique, saturated]
                    out.write('\t'.join(map(str, [c, a+1, b]+values))+'\n')
                    for k, v in zip(['valid_k31_starts','unique_k31_starts','multiple_k31_starts','saturated_frequency_positions'], values):
                        s[k] += v
                seq_summary.append(s)
            assert not f.read(1)
        write_json(OUT / 'reference_k31_opportunity_summary.json', dict(reference_sha256=reference_sha, sequences=seq_summary))
        state('compressing_reference_frequency')
        compressed = raw.with_suffix('.freq8.gz')
        compressed_ck = OUT / 'compression.complete.json'
        if compressed_ck.exists():
            check_records(json.loads(compressed_ck.read_text())['files'])
        else:
            assert not compressed.exists(), 'Existing compressed output requires inspection'
            with raw.open('rb') as f, compressed.open('wb') as dest:
                with gzip.GzipFile(filename='', mode='wb', fileobj=dest, compresslevel=1, mtime=0) as g:
                    for block in iter(lambda: f.read(8*1024*1024), b''):
                        g.write(block)
            # Verify actual decompression, independently of write success.
            h = hashlib.sha256()
            with gzip.open(compressed, 'rb') as f:
                for block in iter(lambda: f.read(8*1024*1024), b''):
                    h.update(block)
            assert h.hexdigest() == sha(raw)
            write_json(compressed_ck, dict(completed_at=stamp(), uncompressed_sha256=h.hexdigest(), files=file_records([compressed])))
        deliver = [compressed, window_path, OUT/'reference_k31_opportunity_summary.json', OUT/'reference_sequence_manifest.json', OUT/'analysis_plan.json', OUT/'toy_v2/validation.json']
        write_json(OUT / 'download_manifest.json', dict(completed_at=stamp(), reference_sha256=reference_sha, files=file_records(deliver), raw_frequency_sha256=sha(raw)))
        state('complete_pending_local_download_and_validation')
    except Exception as e:
        state('failed_preserved_for_checkpoint_inspection', error=repr(e))
        raise
    finally:
        if LOCK.read_text().strip() == str(os.getpid()):
            LOCK.unlink()


if __name__ == '__main__':
    os.umask(0o077)
    os.environ['OMP_NUM_THREADS'] = '8'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    if hasattr(os, 'sched_setaffinity'):
        os.sched_setaffinity(0, set(sorted(os.sched_getaffinity(0))[:8]))
    run()
