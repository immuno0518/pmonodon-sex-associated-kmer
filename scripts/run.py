"""Execute one documented analysis stage with an explicit configuration."""
import argparse, json, os, subprocess, sys
from pathlib import Path

STEPS={
 'prepare':('prepare_run.py',[]),
 'discovery-counts':('run_qc_kmc.py',['--stage','full']),
 'discovery-screen':('screen_discovery_kmers.py',[]),
 'heldout-counts':('validate_heldout_kmers.py',[]),
 'heldout-tests':('analyze_heldout_kmers.py',[]),
 'exact-mapping':('map_heldout_kmers.py',[]),
 'genmap-toy':('test_reference_genmap.py',[]),
 'reference-opportunities':('server_reference_mappability.py',[]),
 'normalized-density':('analyze_reference_mappability.py',[]),
 'existing-depth':('audit_existing_Chr19_depth.py',[]),
 'regional-sensitivity':('analyze_regional_kmer_sensitivity.py',[]),
 'regional-sensitivity-figure':('validate_and_plot_regional_sensitivity.py',[]),
 'figure3-initial':('create_kmer_manuscript_figure.py',[]),
 'figure3-profiles':('revise_kmer_manuscript_figure.py',[]),
 'figure3-discovery-panels':('extract_discovery_figure_panels.py',[]),
 'figure3-all200':('create_all200_kmer_figure.py',[]),
 'figure3-check':('validate_all200_figure.py',[]),
 'regional-vcf':('reviewer_genomics_common.py',[]),
 'local-correlations':('analyze_reviewer_local_ld.py',[]),
 'candidate-sequences':('analyze_reviewer_genes.py',[]),
 'matched-and-genotype-tests':('finalize_reviewer_statistics.py',[]),
 'local-correlation-figure':('plot_reviewer_local_ld.py',[]),
 'secondary-annotations':('secondary_gene_intersections.py',[]),
}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',default='config/local.json')
    ap.add_argument('step',choices=STEPS)
    args,extra=ap.parse_known_args()
    cfg=Path(args.config).resolve()
    if not cfg.is_file():ap.error('Copy config/example.json to config/local.json and configure inputs first.')
    env=os.environ.copy();env['AQUACULTURE_CONFIG']=str(cfg)
    env['PYTHONDONTWRITEBYTECODE']='1'
    if args.step=='discovery-counts':
        c=json.loads(cfg.read_text(encoding='utf-8'));p=Path(c.get('run_root','run'))
        if not p.is_absolute():p=cfg.parent.parent/p
        extra=['--manifest',str(p/'outputs/kmer/discovery_manifest_16.tsv')]+extra
    if args.step=='regional-vcf':
        # The recorded helper masks DP/GQ only when fields are present.
        # Require their definitions before entering the documented quality-mask route.
        import gzip
        c=json.loads(cfg.read_text(encoding='utf-8'))
        source=Path(c['inputs']['variant_source_directory'])
        if not source.is_absolute():source=cfg.parent.parent/source
        for kind in ('snp','indel'):
            p=source/('06.pop_filter_vcf/pop_filter.'+kind+'.vcf.gz')
            with gzip.open(p,'rt') as f:
                definitions=[]
                for line in f:
                    if line.startswith('##FORMAT='):definitions.append(line)
                    if line.startswith('#CHROM'):break
            for field in ('GT','DP','GQ'):
                if not any(('ID='+field+',') in x for x in definitions):
                    raise ValueError('Required VCF FORMAT definition missing: '+field)
    script,base=STEPS[args.step]
    subprocess.run([sys.executable,str(Path(__file__).parent/script),*base,*extra],env=env,check=True)

if __name__=='__main__':main()
