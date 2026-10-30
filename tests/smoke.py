"""Synthetic checks of retained statistical and reference-analysis code.

No study reads, individual identities, original genome or annotation are used.
External KMC, cutadapt and GenMap computations are not exercised here.
"""
import ast, collections, csv, gzip, hashlib, importlib, importlib.metadata
import io, json, math, os, random, runpy, subprocess, sys, unittest
from contextlib import redirect_stdout
from pathlib import Path

PKG=Path(__file__).resolve().parents[1]
WORK=PKG/'run/smoke'
WORK.mkdir(parents=True,exist_ok=True)
CONFIG=WORK/'config/test.json'
CONFIG.parent.mkdir(exist_ok=True)
RUN=WORK/'analysis'
# Paths in this runtime-only configuration are local to the synthetic fixture.
CONFIG.write_text(json.dumps({'run_root':'analysis','reference_sequences':4,'reference_bp':1,
 'inputs':{'reference':'inputs/artificial.fasta','gff':'inputs/artificial.gff3','functional_annotations':'inputs/artificial_annotations.tsv'},
 'annotation_focal_ids':['toy_gene_a']}),encoding='utf-8')
os.environ['AQUACULTURE_CONFIG']=str(CONFIG)
sys.dont_write_bytecode=True
sys.path.insert(0,str(PKG/'scripts'))
import numpy as np
from scipy.stats import beta,fisher_exact
import release_config as cfg
import analyze_heldout_kmers as stats
import map_heldout_kmers as mapper
import reviewer_genomics_common as common
import analyze_reference_mappability as opportunity
import prepare_run

def write_tsv(p,rows,fields=None):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def original_function(file,name,namespace):
    # Load the retained function itself without running unrelated data-dependent stages.
    tree=ast.parse((PKG/'scripts'/file).read_text(encoding='utf-8'))
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
    exec(compile(ast.Module(body=[node],type_ignores=[]),file,'exec'),namespace)
    return namespace[name]

class StatisticalChecks(unittest.TestCase):
    def test_fisher_against_integer_enumeration(self):
        for f,m,nf,nm in [(8,0,8,8),(7,1,8,8),(0,0,8,8),(5,4,9,11),(31,7,92,92),(92,0,92,92),(91,1,92,92)]:
            col=f+m;den=math.comb(nf+nm,nf)
            obs=math.comb(col,f)*math.comb(nf+nm-col,nf-f)
            tail=sum(math.comb(col,x)*math.comb(nf+nm-col,nf-x) for x in range(max(0,nf-(nf+nm-col)),min(nf,col)+1) if math.comb(col,x)*math.comb(nf+nm-col,nf-x)<=obs)
            _,p,corrected=stats.association(f,m,nf,nm,7178221)
            self.assertAlmostEqual(p,tail/den,places=12)
            self.assertEqual(corrected,min(1,p*7178221))
    def test_binomial_interval_endpoints(self):
        lo,hi=stats.proportion_ci(92,92)
        self.assertAlmostEqual(lo,.025**(1/92),places=12);self.assertEqual(hi,1)
        lo,hi=stats.proportion_ci(0,92)
        self.assertEqual(lo,0);self.assertAlmostEqual(hi,1-.025**(1/92),places=12)
    def test_exact_2x3_matches_2x2_when_third_column_empty(self):
        exact23=original_function('finalize_reviewer_statistics.py','exact23',{'math':math})
        self.assertAlmostEqual(exact23([7,1,0],[1,7,0]),stats.association(7,1,8,8,1)[1],places=12)
        self.assertEqual(exact23([2,2,2],[2,2,2]),1)
    def test_sorted_zero_fill_and_unexpected_key(self):
        p=WORK/'counts.tsv';p.write_text('AAA\t2\nTTT\t3\n')
        x=stats.SortedCounts(p)
        self.assertEqual([x.at(k) for k in ['AAA','CCC','TTT']],[2,0,3]);x.finish()
        x=stats.SortedCounts(p)
        try:
            with self.assertRaises(AssertionError):x.at('CCC')
        finally:x.handle.close()
    def test_recorded_genotype_mask(self):
        base=['Chr19','1','.','A','G','.','PASS','.','GT:DP:GQ']
        self.assertEqual(common.genotypes(base+['0/1:5:20','1|1:9:30','0/0:4:99','0/1:9:19','./.:.:.','1/2:9:30']),[1,2,-1,-1,-1,-1])
    def test_single_variant_consequence_and_splice(self):
        from Bio.Seq import Seq
        reference=lambda c,a,b:'ATGAAATAA'[a:b]
        fn=original_function('analyze_reviewer_genes.py','consequence',{'reference':reference,'Seq':Seq})
        g={'strand':'+','map':{i:i-1 for i in range(1,10)},'cds':'ATGAAATAA','cds_exons':[(1,9,0)]}
        self.assertEqual(fn(['Chr19','4','.','A','G'],g)[0],'missense')
        self.assertEqual(fn(['Chr19','4','.','A','AT'],g)[0],'predicted_frameshift')
        self.assertEqual(fn(['Chr19','6','.','A','G'],g)[0],'synonymous')
        g['map']={i:i-1 for i in range(1,4)};g['cds_exons']=[(1,3,0),(10,12,0)]
        self.assertEqual(fn(['Chr19','4','.','A','G'],g)[0],'predicted_splice_donor')
    def test_source_self_checks(self):
        stats.self_test();mapper.self_test()

class ReferenceAndAnnotationChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=random.Random(47)
        keys=[]
        while len(keys)<4:
            k=''.join(rng.choice('ACGT') for _ in range(31));k=min(k,k.translate(mapper.RC)[::-1])
            if k not in keys:keys.append(k)
        cls.keys=keys
        a,b,c,d=keys
        cls.seqs={'Chr19':'NN'+a+'NN'+b.translate(mapper.RC)[::-1]+'NN',
                  'Chr34':c+'N'+c,'Contig6':'NN'+b+'NN','short':'ACGT'}
        cls.seqs['Contig6']='NN'+b+'NN'
        cls.keys=[a,d,c,b]
        cls.unmapped='T'*31
        cls.seqs['Chr19']='NN'+a+'NN'+d.translate(mapper.RC)[::-1]+'NN'
        cfg.CONFIG['reference_bp']=sum(map(len,cls.seqs.values()))
        cls.config=json.loads(CONFIG.read_text());cls.config['reference_bp']=cfg.CONFIG['reference_bp']
        CONFIG.write_text(json.dumps(cls.config))
        for p in [RUN/'inputs',mapper.WORK,mapper.RESULT,opportunity.DOWNLOAD,opportunity.BASE,WORK/'inputs']:
            p.mkdir(parents=True,exist_ok=True)
        ref=WORK/'inputs/artificial.fasta'
        ref.write_text(''.join('>'+name+'\n'+s+'\n' for name,s in cls.seqs.items()))
        entries=[dict(original_id=n,length_bp=len(s),sequence_sha256=hashlib.sha256(s.encode()).hexdigest()) for n,s in cls.seqs.items()]
        write_tsv(RUN/'inputs/reference_sequence_manifest.tsv',entries)
        cls.entries=entries
        cls.queries=[dict(query_id=f'toy_{i}',candidate_group=g,kmer=k,heldout_F_present=92 if g=='F' else 0,heldout_M_present=92 if g=='M' else 0,bonferroni_p='1e-20',perfect_heldout_separation='1') for i,(g,k) in enumerate([('F',a),('M',d),('F',c),('M',cls.unmapped)])]
        qp=mapper.WORK/'queries_000.tsv';write_tsv(qp,cls.queries,mapper.QUERY_COLUMNS)
        cls.item=dict(file=qp.name,count=4,sha256=mapper.digest(qp))
        cls.manifest={'chunks':[cls.item],'total_queries':4}
        # Frequency8 fixtures are independently enumerated Python values, not GenMap outputs.
        freq=collections.Counter()
        for s in cls.seqs.values():
            for i in range(max(0,len(s)-30)):
                k=s[i:i+31]
                if set(k)<=set('ACGT'):freq[min(k,k.translate(mapper.RC)[::-1])]+=1
        raw=bytes(min(255,freq[min(s[i:i+31],s[i:i+31].translate(mapper.RC)[::-1])]) if len(s[i:i+31])==31 and set(s[i:i+31])<=set('ACGT') else 0 for s in cls.seqs.values() for i in range(len(s)))
        cls.raw=raw
        raw_sha=hashlib.sha256(raw).hexdigest()
        with gzip.open(opportunity.DOWNLOAD/'data01_k31_e0.freq8.gz','wb') as f:f.write(raw)
        records=[dict(sequence=e['original_id'],length_bp=e['length_bp'],sequence_sha256=e['sequence_sha256']) for e in entries]
        (opportunity.DOWNLOAD/'reference_sequence_manifest.json').write_text(json.dumps({'sequences':records}))
        windows=[];summary=[];offset=0
        for e in entries:
            data=raw[offset:offset+e['length_bp']];offset+=e['length_bp']
            valid=len(data)-data.count(0);unique=data.count(1)
            row=dict(sequence=e['original_id'],start_1based=1,end_1based=e['length_bp'],valid_k31_starts=valid,unique_k31_starts=unique,multiple_k31_starts=valid-unique,saturated_frequency_positions=data.count(255))
            windows.append(row);summary.append(dict(sequence=e['original_id'],valid_k31_starts=valid,unique_k31_starts=unique))
        write_tsv(opportunity.DOWNLOAD/'reference_k31_opportunity_100kb.tsv',windows)
        (opportunity.DOWNLOAD/'reference_k31_opportunity_summary.json').write_text(json.dumps({'sequences':summary}))
        paths=[p for p in opportunity.DOWNLOAD.iterdir() if p.is_file() and p.name!='download_manifest.json']
        (opportunity.DOWNLOAD/'download_manifest.json').write_text(json.dumps({'files':[dict(path=p.name,bytes=p.stat().st_size,sha256=mapper.digest(p)) for p in paths],'raw_frequency_sha256':raw_sha}))
    def test_01_full_exact_mapping_and_zero_one_based_edges(self):
        # Use the actual chunk mapper, writer and merged-result code.
        mapper.map_chunk(self.item,0,1,mapper.reference_identity())
        result=mapper.merge_results(self.manifest)
        self.assertEqual(result['categories'],{'F':{'unique':1,'multiple':1},'M':{'unique':1,'unmapped':1}})
        with gzip.open(mapper.RESULT/'all_replicated_candidates_reference_mapping.tsv.gz','rt') as f:r=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual([int(x['exact_reference_hits']) for x in r],[1,1,2,0])
        self.assertEqual((r[0]['unique_start_1based'],r[0]['unique_end_1based'],r[0]['unique_strand']),('3','33','+'))
        self.assertIn('Chr19:36-66:-',r[1]['first_3_hits'])
    def test_02_unique_opportunity_consumer(self):
        # Supply a minimal consistent mapping fixture with both groups represented.
        # All rows are independently checked against the actual reference sequence.
        unique=[]
        observed={}
        for name,seq in self.seqs.items():
            for i in range(max(0,len(seq)-30)):
                k=seq[i:i+31]
                if set(k)<=set('ACGT'):
                    key=min(k,k.translate(mapper.RC)[::-1]);observed.setdefault(key,[]).append((name,i+1,'+' if key==k else '-'))
        singles=[(k,h[0]) for k,h in observed.items() if len(h)==1]
        self.assertGreaterEqual(len(singles),2)
        for i,(k,(name,start,strand)) in enumerate(singles[:2]):
            g='F' if i==0 else 'M'
            unique.append(dict(query_id='norm_'+str(i),candidate_group=g,kmer=k,heldout_F_present=92,heldout_M_present=0,bonferroni_p='1e-20',perfect_heldout_separation='1',exact_reference_hits=1,unique_chr=name,unique_start_1based=start,unique_end_1based=start+30,unique_strand=strand,first_3_hits=f'{name}:{start}-{start+30}:{strand}',mononucleotide_entropy_bits='2'))
        with gzip.open(mapper.RESULT/'all_replicated_candidates_reference_mapping.tsv.gz','wt',newline='') as f:
            w=csv.DictWriter(f,fieldnames=mapper.MAP_COLUMNS,delimiter='\t');w.writeheader();w.writerows(unique)
        (mapper.RESULT/'mapping_summary.json').write_text(json.dumps({'categories':{'F':{'unique':1},'M':{'unique':1}}}))
        opportunity.main()
        result=json.loads((opportunity.BASE/'opportunity_normalized_regional_assessment.json').read_text())
        self.assertEqual(result['reference_unique_k31_starts'],self.raw.count(1))
        self.assertEqual(result['total_associated_unique'],{'F':1,'M':1})
        self.assertEqual(result['target']['unique_k31_starts'],0)
        self.assertIsNone(result['target']['F_associated_per_million_unique_starts'])
    def test_03_annotation_intersection_stage(self):
        rows=[dict(candidate_group=g,chromosome=c,start_1based=1,end_1based=40,kmer_count=n) for g,c,n in [('F','Chr34',5),('F','Chr34',2),('M','Contig6',7)]]
        write_tsv(mapper.RESULT/'unique_mapping_clusters_gap1kb.tsv',rows)
        (WORK/'inputs/artificial.gff3').write_text('##gff-version 3\nChr34\ttoy\tgene\t10\t60\t.\t+\t.\tID=toy_gene_a\nContig6\ttoy\tgene\t1\t31\t.\t+\t.\tID=toy_gene_b\n')
        fields=['Gene_id','Chrom','Gene_start','Gene_end','Strand','Query_length','NR_Accession','NR_Description','SwissProt_Gene_Name','Pfam_IDs','Pfam_Descriptions']
        write_tsv(WORK/'inputs/artificial_annotations.tsv',[{k:('toy_gene_a' if k=='Gene_id' else 'synthetic') for k in fields}],fields)
        runpy.run_path(str(PKG/'scripts/secondary_gene_intersections.py'),run_name='__main__')
        result=json.loads((RUN/'outputs/annotation/secondary_cluster_gene_annotation.json').read_text())
        self.assertEqual(sorted(int(x['cluster_kmers']) for x in result['overlapping_genes']),[5,7])
        self.assertEqual(next(x for x in result['overlapping_genes'] if x['gene_id']=='toy_gene_b')['NR_Description'],'not available')

class LocalCorrelationChecks(unittest.TestCase):
    def test_synthetic_200_individual_correlation_stage(self):
        ids=['toy%03d'%i for i in range(200)]
        write_tsv(RUN/'outputs/kmer/sample_manifest_200.tsv',[dict(sample_id=s,sex='F' if i<100 else 'M') for i,s in enumerate(ids)])
        common.OUT.mkdir(parents=True,exist_ok=True)
        p=common.OUT/'Chr19_2500001_3100000.snp.vcf.gz'
        patterns=[[i%3 for i in range(200)],[(i//3)%3 for i in range(200)],[(i//7)%3 for i in range(200)],
                  [1 if i<100 else 0 for i in range(200)],[1]*200,[0 if i<100 else i%3 for i in range(200)]]
        with gzip.open(p,'wt') as f:
            f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n')
            for j,values in enumerate(patterns):
                gt=['0/0','0/1','1/1']
                f.write('Chr19\t%d\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t'%(2500001+j*10000)+'\t'.join(gt[v]+':8:30' for v in values)+'\n')
        runpy.run_path(str(PKG/'scripts/analyze_reviewer_local_ld.py'),run_name='__main__')
        out=common.OUT/'local_LD';summary=json.loads((out/'summary.json').read_text())
        self.assertGreater(summary['independent_direct_pair_checks'],0)
        self.assertLess(summary['max_absolute_error'],1e-10)
        self.assertGreater(summary['results']['5']['groups']['female100']['constant_genotype_sites_excluded'],0)
        for group in ('all200','female100','male100'):
            d=np.load(out/(group+'_DP5_pairwise_r2.npz'))
            self.assertEqual(d['r2'].shape,(len(d['positions']),len(d['positions'])))
            finite=d['r2'][np.isfinite(d['r2'])]
            self.assertTrue(np.all((finite>=0)&(finite<=1)))

class PackagingChecks(unittest.TestCase):
    def test_stage_dispatcher_creates_separate_run_directories(self):
        p=subprocess.run([sys.executable,str(PKG/'scripts/run.py'),'--config',str(CONFIG),'prepare'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertTrue((RUN/'work/kmer/discovery').is_dir())
    def test_all_python_parse(self):
        for p in list((PKG/'scripts').glob('*.py'))+list((PKG/'tests').glob('*.py')):
            ast.parse(p.read_text(encoding='utf-8'),filename=p.name)
    def test_source_hash_map_matches_files(self):
        for r in json.loads((PKG/'provenance/script_sources.json').read_text()):
            self.assertEqual(hashlib.sha256((PKG/r['released_file']).read_bytes()).hexdigest(),r['released_sha256'])
    def test_input_manifest_rejects_bad_sex(self):
        p=WORK/'bad.tsv';write_tsv(p,[dict(sample_id='artificial',sex='unknown',R1='r1',R2='r2',sequencing_batch='Not_recorded')])
        with self.assertRaises(ValueError):prepare_run.read_samples(p)
    def test_aggregate_source_denominators(self):
        p=PKG/'source_data/Figure3_panel_A_chromosome_density.tsv'
        with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(sum(int(x['unique_k31_starts']) for x in rows),1139681763)
        self.assertEqual([sum(int(x[g+'_unique']) for x in rows) for g in 'FM'],[180082,184371])
        for x in rows:
            for g in 'FM':self.assertAlmostEqual(float(x[g+'_per_1000_unique_starts']),1000*int(x[g+'_unique'])/int(x['unique_k31_starts']),places=9)

def main():
    suite=unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    # Console paths are local test execution details; the archived report contains names only.
    output=io.StringIO();capture=io.StringIO()
    with redirect_stdout(capture):result=unittest.TextTestRunner(stream=output,verbosity=2).run(suite)
    env={'python':sys.version.split()[0]}
    for p in ['numpy','scipy','biopython','pyahocorasick','psutil','matplotlib']:
        try:env[p]=importlib.metadata.version(p)
        except importlib.metadata.PackageNotFoundError:env[p]='not installed'
    report={'status':'passed' if result.wasSuccessful() else 'failed','tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'environment':env,
     'fixture':'Artificial sequences, generated canonical frequencies, mock annotation and simple genotype records; no study data',
     'external_tools_executed':[],'study_full_rerun':False,'genmap_frequency_fixture':'Independent Python enumeration, not a GenMap execution',
     'results':output.getvalue().splitlines()}
    (PKG/'validation/test_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(output.getvalue());print(json.dumps({k:v for k,v in report.items() if k!='results'},indent=2))
    return 0 if result.wasSuccessful() else 1

if __name__=='__main__':sys.exit(main())
