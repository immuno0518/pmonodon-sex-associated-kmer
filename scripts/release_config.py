"""Configuration for the portable copies of the recorded analysis scripts."""
import csv, hashlib, json, os, shutil
from pathlib import Path

CONFIG_FILE = Path(os.environ.get('AQUACULTURE_CONFIG', 'config/local.json')).resolve()
CONFIG = json.loads(CONFIG_FILE.read_text(encoding='utf-8')) if CONFIG_FILE.exists() else {}
PROJECT = CONFIG_FILE.parent.parent
def resolve_input(value):
    p=Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (PROJECT/p).resolve()
RUN_ROOT = resolve_input(CONFIG.get('run_root','run'))
def input_path(key):
    return resolve_input(CONFIG.get('inputs',{}).get(key,'inputs/'+key))
def binary(key,default):
    v=CONFIG.get('tools',{}).get(key,default)
    if '/' in v or chr(92) in v:return resolve_input(v)
    return Path(shutil.which(v) or v)
KMC_BIN=binary('kmc','kmc')
KMC_TOOLS_BIN=binary('kmc_tools','kmc_tools')
GENMAP_BIN=binary('genmap','genmap')
def verify_binary(path,key):
    expected=CONFIG.get('binary_checksums',{}).get(key)
    if expected:
        h=hashlib.sha256(Path(path).read_bytes()).hexdigest()
        if h!=expected:raise ValueError('Executable checksum mismatch: '+key)
def expected_sequences():return CONFIG.get('reference_sequences',93)
def expected_reference_bp():return CONFIG.get('reference_bp',2111532351)
def identity_names():
    with (RUN_ROOT/'inputs/reference_sequence_manifest.tsv').open(encoding='utf-8') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    return {r['original_id']:dict(r,new_id=r['original_id']) for r in rows}
