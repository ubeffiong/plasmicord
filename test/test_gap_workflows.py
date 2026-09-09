"""Regression coverage for annotation, native boundaries, quality and calibration."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from plasmid_annotation import annotate, parse_gff, parse_amrfinder, to_plasbench
from plasmid_annotation.engine import database_identity
from python.importers import native_candidates, import_inputs
from python.quality import assess, load_evidence, graph_closure
from python.contracts import write_tsv, read_tsv
from python.calibration import calibrate
from python.cli import aggregate, run
from python.external_typing import load_external_typing
from python.containment import detect_containment
from python.population_summary import summarize
from python.multilayer_network import cluster_relation, multilayer_edges


class GapWorkflows(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def gff(self, name='contig1'):
        return f'##gff-version 3\n{name}\tBakta\tCDS\t10\t90\t.\t-\t0\tID=cds1;gene=merA;product=mercuric%20reductase\n##FASTA\n>{name}\nACGT\n'

    def test_gff_coordinates_category_and_plasbench_export(self):
        f=self.p/'features.gff';f.write_text(self.gff())
        row=parse_gff(f,'bakta','1.12','bakta','6')[0]
        self.assertEqual(row['functional_category'],'metal_resistance')
        self.assertEqual(row['start'],10)
        self.assertEqual(to_plasbench([row])[0]['start'],9)
        self.assertEqual(to_plasbench([row])[0]['end'],90)
        self.assertNotIn('amr_gene',row)

    def test_partial_amrfinder_hits_are_not_headlines(self):
        f=self.p/'amr.tsv'
        fields=['Contig id','Start','Stop','Strand','Gene symbol','Sequence name','Element type','Element subtype','Method','% Identity to reference sequence','% Coverage of reference sequence','Class']
        rows=[dict(zip(fields,['contig1',1,90,'+','blaTEST','test gene','AMR','AMR','EXACT',100,100,'BETA-LACTAM'])),
              dict(zip(fields,['contig1',100,190,'-','blaPART','test partial','AMR','AMR','PARTIAL',100,100,'BETA-LACTAM']))]
        write_tsv(f,fields,rows)
        parsed=parse_amrfinder(f,'4.2','2026')
        self.assertEqual(parsed[0]['hit_class'],'strict')
        self.assertEqual(parsed[1]['hit_class'],'partial')
        self.assertEqual(parsed[1]['annotation_confidence'],'low')

    def test_cache_binds_sequences_to_new_ids_and_invalidates_database_changes(self):
        db=self.p/'db';db.mkdir();(db/'VERSION').write_text('v1')
        config={'genes':dict(engine='bakta',database=str(db),database_version='1')}
        calls=[]
        def fake_run(cmd,**kwargs):
            if '--version' in cmd:return subprocess.CompletedProcess(cmd,0,stdout='bakta 1.12',stderr='')
            calls.append(cmd)
            dest=Path(cmd[cmd.index('--output')+1]);dest.mkdir()
            (dest/'annotation.gff3').write_text(self.gff())
            return subprocess.CompletedProcess(cmd,0)
        with patch('plasmid_annotation.engine.shutil.which',return_value='/mock/bakta'),patch('plasmid_annotation.engine.subprocess.run',side_effect=fake_run):
            first,_,_=annotate([('original','ACGT'*100)],config,self.p/'one',self.p/'cache')
            second,status,_=annotate([('renamed','ACGT'*100)],config,self.p/'two',self.p/'cache')
            self.assertEqual(second[0]['source_sequence_id'],'renamed')
            self.assertEqual(status[0]['cache'],'reused')
            self.assertEqual(len(calls),1)
            (db/'VERSION').write_text('changed content with same declared version')
            annotate([('renamed','ACGT'*100)],config,self.p/'three',self.p/'cache')
            self.assertEqual(len(calls),2)

    def test_modern_amrfinder_schema_and_plus_scope(self):
        fields=['Contig id','Start','Stop','Strand','Element symbol','Element name','Type','Subtype','Method','% Identity to reference','% Coverage of reference','Scope']
        values=['c',1,99,'+','blaTEST','test','AMR','AMR','EXACTX',100,100,'core']
        rows=[dict(zip(fields,values))]
        rows.append(dict(rows[0], Scope='plus'))
        rows.append(dict(rows[0], **{'% Identity to reference':'NA','Method':'HMM'}))
        f=self.p/'modern.tsv';write_tsv(f,fields,rows)
        parsed=parse_amrfinder(f,'4.2.7','2026-08-07.1')
        self.assertEqual(parsed[0]['hit_class'],'strict')
        self.assertEqual(parsed[1]['hit_class'],'supporting')
        self.assertEqual(parsed[2]['identity'],'')

    def test_calibration_rejects_mismatched_engine_provenance(self):
        args,_=self.calibration_fixture()
        (self.p/'run_provenance.json').write_text(json.dumps(dict(status='complete',distance_engine='kmer')))
        with self.assertRaisesRegex(ValueError,'engine must match'):
            calibrate(args)
        (self.p/'run_provenance.json').write_text(json.dumps(dict(status='complete',distance_engine='mash',
            output_sha256={args.index.name:'0'*64})))
        with self.assertRaisesRegex(ValueError,'index checksum'):
            calibrate(args)

    def test_annotation_failure_is_not_a_negative_result(self):
        db=self.p/'db';db.mkdir();(db/'v').write_text('1')
        config={'genes':dict(engine='bakta',database=str(db),database_version='1')}
        with patch('plasmid_annotation.engine.shutil.which',return_value=None):
            with self.assertRaisesRegex(ValueError,'unavailable'):
                annotate([('c','ACGT')],config,self.p/'out')

    def test_flye_requires_selection_and_preserves_circularity_level(self):
        folder=self.p/'flye';folder.mkdir()
        (folder/'assembly.fasta').write_text('>contig_1\nACGTACGT\n>contig_2\nTGCATGCA\n')
        (folder/'assembly_info.txt').write_text('#seq_name\tlength\tcov.\tcirc.\ncontig_1\t8\t30\tY\ncontig_2\t8\t20\tN\n')
        with self.assertRaisesRegex(ValueError,'Explicit'):
            list(native_candidates('flye',folder))
        rows=list(native_candidates('flye',folder,'contig_1'))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0][2]['circularity_status'],'tool_reported')

    def test_tadrep_preserves_alignment_provenance_and_leaves_circularity_unresolved(self):
        folder=self.p/'tadrep';folder.mkdir()
        (folder/'sample1-refA-pseudo.fna').write_text('>refA\nACGTACGTACGT\n')
        header='plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length'
        (folder/'sample1-summary.tsv').write_text(header+'\n'
            'refA\tcontig_1\t1\t8\t8\t95.0\t99.0\t8\t+\t1\t8\t12\n'
            'refA\tcontig_2\t1\t4\t4\t90.0\t97.5\t4\t+\t9\t12\t12\n')
        rows=list(native_candidates('tadrep',folder))
        self.assertEqual(len(rows),1)
        name,records,extra,evidence=rows[0]
        self.assertEqual(name,'refA')
        self.assertEqual(extra['circularity_status'],'unresolved')
        self.assertEqual(extra['reported_coverage'],'90.0')
        self.assertEqual(extra['reported_identity'],'97.5')
        self.assertEqual(extra['alignment_length_bp'],'12')
        (folder/'sample1-summary.tsv').unlink()
        with self.assertRaisesRegex(ValueError,'summary.tsv'):
            list(native_candidates('tadrep',folder))

    def test_tadrep_rejects_pseudo_file_with_no_matching_summary_row(self):
        folder=self.p/'tadrep_mismatch';folder.mkdir()
        (folder/'sample1-refB-pseudo.fna').write_text('>refB\nACGTACGTACGT\n')
        header='plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length'
        (folder/'sample1-summary.tsv').write_text(header+'\nrefA\tcontig_1\t1\t8\t8\t95.0\t99.0\t8\t+\t1\t8\t12\n')
        with self.assertRaisesRegex(ValueError,'no matching rows'):
            list(native_candidates('tadrep',folder))

    def test_completion_denominator_tracks_engine_not_product_category(self):
        index=[dict(plasmid_id='p',isolate_id='i')]
        features=[dict(plasmid_id='p',isolate_id='i',plasmid_unit='PU',feature_id='f',
                       functional_category='metal_resistance',gene_symbol='merA',annotation_engine='prokka')]
        statuses=[dict(plasmid_id='p',stage='genes',status='complete')]
        row=aggregate(index,features,[dict(isolate_id='i')],{'p':'PU'},statuses)[0]
        self.assertEqual(row['n_evaluated_plasmids'],1)
        path=self.p/'single.fasta';path.write_text('>one\nACGT\n')
        self.assertEqual(len(list(native_candidates('generic-fasta',path))),1)

    def test_mob_contig_boundaries_and_chromosome_exclusion(self):
        (self.p/'plasmid_X.fasta').write_text('>a\nACGT\n>b\nTGCA\n')
        (self.p/'contig_report.txt').write_text('contig_id\tmolecule_type\na\tplasmid\nb\tplasmid\n')
        rows=list(native_candidates('mob-recon',self.p))
        self.assertEqual(len(rows[0][1]),2)
        (self.p/'contig_report.txt').write_text('contig_id\tmolecule_type\na\tplasmid\nb\tchromosome\n')
        with self.assertRaisesRegex(ValueError,'plasmid assignments'):
            list(native_candidates('mob-recon',self.p))

    def test_plasbench_adapter_retains_selection_and_requires_candidate_boundaries(self):
        (self.p/'selection_report.json').write_text(json.dumps(dict(selected_tool='mob_recon',selection_type='truth_set_best_candidate')))
        (self.p/'candidate.plasmid.fasta').write_text('>a\nACGT\n>b\nTGCA\n')
        with self.assertRaisesRegex(ValueError,'ambiguous plasmid boundaries'):
            list(native_candidates('plasbench',self.p))
        args=argparse.Namespace(adapter='plasbench',input=self.p,isolate_id='iso',samples=None,record_ids='a',group_contigs=False,tool_version='1',out=self.p/'imported')
        import_inputs(args)
        self.assertIn('mob_recon',(args.out/'manifest.tsv').read_text())
        self.assertTrue((args.out/'import_provenance.json').exists())

    def test_unicycler_and_generic_do_not_assume_all_contigs_are_plasmids(self):
        (self.p/'assembly.fasta').write_text('>1 length=8 circular=true\nACGTACGT\n>2\nTGCATGCA\n')
        row=list(native_candidates('unicycler',self.p,'1'))[0]
        self.assertEqual(row[2]['circularity_status'],'tool_reported')
        with self.assertRaisesRegex(ValueError,'Explicit'):
            list(native_candidates('generic-fasta',self.p/'assembly.fasta'))

    def quality_index(self,sequence='ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC'):
        return [dict(plasmid_id='p1',sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),quality_status='uncertain',quality_warnings='',circularity_status='unresolved')],{'p1':[('c1',sequence)]}

    def test_quality_does_not_reject_novel_untyped_plasmids(self):
        index,seq=self.quality_index()
        row=assess(index,seq)[0]
        self.assertEqual(row['quality_status'],'uncertain')
        self.assertEqual(row['classification'],'uncertain')

    def test_high_quality_needs_sequence_bound_graph_and_independent_evidence(self):
        index,seq=self.quality_index()
        graph=self.p/'graph.gfa';graph.write_text(f"S\tc1\t{seq['p1'][0][1]}\nL\tc1\t+\tc1\t+\t0M\n")
        index[0]['assembly_graph_path']=str(graph)
        evidence={'p1':dict(classification='plasmid',read_breadth='1',mean_depth='30',chromosome_fraction='0',evidence_source='reference and reads')}
        self.assertEqual(assess(index,seq,evidence=evidence)[0]['quality_status'],'high_confidence')
        graph.write_text('S\tc1\tAAAA\nL\tc1\t+\tc1\t+\t0M\n')
        self.assertEqual(graph_closure(graph,seq['p1']),'unresolved')

    def test_chromosome_evidence_rejects_and_checksum_mismatch_fails(self):
        index,seq=self.quality_index()
        self.assertEqual(assess(index,seq,evidence={'p1':dict(classification='chromosome')})[0]['quality_status'],'rejected')
        path=self.p/'evidence.tsv'
        write_tsv(path,['plasmid_id','sequence_sha256','evidence_source'],[dict(plasmid_id='p1',sequence_sha256='wrong',evidence_source='lab')])
        with self.assertRaisesRegex(ValueError,'checksum'):
            load_evidence(path,index)

    def calibration_fixture(self):
        ids=list('ABCDEFGH')
        matrix=[[0 if i==j else .5 for j in range(8)] for i in range(8)]
        matrix[0][1]=matrix[1][0]=.005
        matrix[4][5]=matrix[5][4]=.02
        f=self.p/'matrix.tsv';f.write_text('item\t'+'\t'.join(ids)+'\n'+''.join(a+'\t'+'\t'.join(map(str,row))+'\n' for a,row in zip(ids,matrix)))
        index=self.p/'index.tsv'
        write_tsv(index,['plasmid_id','isolate_id','sequence_sha256'],[dict(plasmid_id=i,isolate_id=i,sequence_sha256=hashlib.sha256(i.encode()).hexdigest()) for i in ids])
        labels=self.p/'labels.tsv'
        rows=[dict(plasmid_a=a,plasmid_b=b,label=label,split=split,group_id=split,evidence_source='independent laboratory evidence',evidence_type='experimental') for a,b,label,split in [('A','B','positive','train'),('C','D','negative','train'),('E','F','positive','validation'),('G','H','negative','validation')]]
        write_tsv(labels,list(rows[0]),rows)
        return argparse.Namespace(matrix=f,index=index,labels=labels,cohort_kind='real',engine='mash',thresholds='0.01,0.03',linkage='complete',min_groups=1,out=self.p/'calibration'),rows

    def test_calibration_selects_training_threshold_without_peeking(self):
        args,_=self.calibration_fixture();calibrate(args)
        result=json.loads((args.out/'calibration.json').read_text())
        self.assertEqual(result['selected_threshold'],.01)
        self.assertIn('validation\t0.01\t0\t0\t1\t1',(args.out/'calibration_metrics.tsv').read_text())

    def test_calibration_abstains_when_group_support_insufficient(self):
        args,_=self.calibration_fixture();args.min_groups=3;calibrate(args)
        self.assertIsNone(json.loads((args.out/'calibration.json').read_text())['selected_threshold'])

    def test_calibration_rejects_group_or_sequence_leakage(self):
        args,rows=self.calibration_fixture()
        rows[-1]['group_id']='train';write_tsv(args.labels,list(rows[0]),rows)
        with self.assertRaisesRegex(ValueError,'group spans'):
            calibrate(args)

    def test_external_typing_rejects_checksum_mismatch_and_unknown_id(self):
        index,_=self.quality_index()
        fields=['plasmid_id','sequence_sha256','evidence_source','external_tool']
        path=self.p/'typing.tsv'
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256='wrong',evidence_source='lab',external_tool='mob_suite')])
        with self.assertRaisesRegex(ValueError,'checksum'):
            load_external_typing(path,index)
        write_tsv(path,fields,[dict(plasmid_id='unknown',sequence_sha256='x',evidence_source='lab',external_tool='mob_suite')])
        with self.assertRaisesRegex(ValueError,'checksum'):
            load_external_typing(path,index)

    def test_external_typing_rejects_blank_external_tool(self):
        index,_=self.quality_index()
        digest=index[0]['sequence_sha256']
        fields=['plasmid_id','sequence_sha256','evidence_source','external_tool']
        path=self.p/'typing.tsv'
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='')])
        with self.assertRaisesRegex(ValueError,'external_tool'):
            load_external_typing(path,index)

    def test_external_typing_round_trips_opaque_fields_and_validates_confidence(self):
        index,_=self.quality_index()
        digest=index[0]['sequence_sha256']
        fields=['plasmid_id','sequence_sha256','evidence_source','external_tool','mob_primary_cluster_id','ptu_assignment','ptu_confidence','associated_pmids']
        path=self.p/'typing.tsv'
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='mob_suite',
                                    mob_primary_cluster_id='AA123',ptu_assignment='PTU-FE',ptu_confidence='high',associated_pmids='12345;67890')])
        result=load_external_typing(path,index)
        self.assertEqual(result['p1']['mob_primary_cluster_id'],'AA123')
        self.assertEqual(result['p1']['ptu_assignment'],'PTU-FE')
        self.assertEqual(result['p1']['associated_pmids'],'12345;67890')
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='mob_suite',
                                    mob_primary_cluster_id='',ptu_assignment='',ptu_confidence='not-a-number',associated_pmids='')])
        with self.assertRaisesRegex(ValueError,'ptu_confidence'):
            load_external_typing(path,index)
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='mob_suite',
                                    mob_primary_cluster_id='',ptu_assignment='',ptu_confidence='',associated_pmids='not-a-pmid')])
        with self.assertRaisesRegex(ValueError,'associated_pmids'):
            load_external_typing(path,index)

    def test_external_typing_validates_predicted_transmissibility_fields(self):
        index,_=self.quality_index()
        digest=index[0]['sequence_sha256']
        fields=['plasmid_id','sequence_sha256','evidence_source','external_tool','predicted_transmissibility_score','predicted_transmissibility_call']
        path=self.p/'typing.tsv'
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='plastrans',
                                    predicted_transmissibility_score='0.87',predicted_transmissibility_call='conjugative')])
        result=load_external_typing(path,index)
        self.assertEqual(result['p1']['predicted_transmissibility_score'],'0.87')
        self.assertEqual(result['p1']['predicted_transmissibility_call'],'conjugative')
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='plastrans',
                                    predicted_transmissibility_score='1.5',predicted_transmissibility_call='')])
        with self.assertRaisesRegex(ValueError,'predicted_transmissibility_score'):
            load_external_typing(path,index)
        write_tsv(path,fields,[dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='plastrans',
                                    predicted_transmissibility_score='',predicted_transmissibility_call='definitely-transferable')])
        with self.assertRaisesRegex(ValueError,'predicted_transmissibility_call'):
            load_external_typing(path,index)

    def test_containment_heuristic_skips_equal_length_and_size_disparity(self):
        accepted=[dict(plasmid_id='small',isolate_id='i1',length=950),
                  dict(plasmid_id='large',isolate_id='i2',length=1000),
                  dict(plasmid_id='equal_a',isolate_id='i3',length=1000),
                  dict(plasmid_id='equal_b',isolate_id='i4',length=1000),
                  dict(plasmid_id='tiny',isolate_id='i5',length=100)]
        ids=[r['plasmid_id'] for r in accepted]
        pos={pid:i for i,pid in enumerate(ids)}
        n=len(ids);distances=[[0.0]*n for _ in range(n)]
        def setd(a,b,v): distances[pos[a]][pos[b]]=distances[pos[b]][pos[a]]=v
        setd('small','large',0.02);setd('equal_a','equal_b',0.01);setd('tiny','large',0.01)
        rows=detect_containment(accepted,ids,pos,distances)
        pairs={(r['small_plasmid_id'],r['large_plasmid_id']) for r in rows}
        self.assertIn(('small','large'),pairs)
        self.assertNotIn(('equal_a','equal_b'),pairs)
        self.assertNotIn(('tiny','large'),pairs)
        for r in rows:
            self.assertEqual(r['heuristic_flag'],'heuristic_length_similarity')
            self.assertIn('NOT alignment-confirmed',r['interpretation'])

    def test_containment_rejects_invalid_ratio_bounds(self):
        with self.assertRaisesRegex(ValueError,'containment-min-ratio'):
            detect_containment([],[],{},[],min_length_ratio=0.9,max_length_ratio=0.5)

    def run_fixture(self,out_name,external_typing=None,features=None):
        meta=self.p/'meta.tsv'
        write_tsv(meta,['isolate_id'],[dict(isolate_id='i1'),dict(isolate_id='i2')])
        seq='ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC'
        fasta1=self.p/'p1.fa';fasta1.write_text(f'>p1\n{seq}\n')
        fasta2=self.p/'p2.fa';fasta2.write_text(f'>p2\n{seq}\n')
        manifest=self.p/f'{out_name}_manifest.tsv'
        write_tsv(manifest,['isolate_id','plasmid_id','fasta_path'],
                  [dict(isolate_id='i1',plasmid_id='p1',fasta_path=str(fasta1)),
                   dict(isolate_id='i2',plasmid_id='p2',fasta_path=str(fasta2))])
        args=argparse.Namespace(manifest=manifest,metadata=meta,features=features,mode='precomputed',engine='kmer',
            out=self.p/out_name,threshold=.05,k=3,min_length=3,sketch_size=100,linkage='complete',
            annotation_config=None,external_typing=external_typing)
        run(args)
        return args.out

    def test_cluster_relation_labels_same_cross_and_unknown(self):
        meta_by_iso={'i1':dict(chromosomal_cluster='CC1'),'i2':dict(chromosomal_cluster='CC1'),
                     'i3':dict(chromosomal_cluster='CC2'),'i4':dict(chromosomal_cluster='')}
        self.assertEqual(cluster_relation(meta_by_iso,'i1','i2'),'same_cluster')
        self.assertEqual(cluster_relation(meta_by_iso,'i1','i3'),'cross_cluster')
        self.assertEqual(cluster_relation(meta_by_iso,'i1','i4'),'unknown_cluster')
        self.assertEqual(cluster_relation(meta_by_iso,'i3','i4'),'unknown_cluster')

    def test_multilayer_edges_does_not_aggregate_across_shared_units(self):
        meta_by_iso={'i1':dict(chromosomal_cluster='CC1'),'i2':dict(chromosomal_cluster='CC1')}
        edge_details=[dict(source='i1',target='i2',plasmid_unit='PU_0001',minimum_distance=0.0,threshold_margin=.05),
                      dict(source='i1',target='i2',plasmid_unit='PU_0002',minimum_distance=.01,threshold_margin=.04)]
        layered=multilayer_edges(edge_details,meta_by_iso)
        self.assertEqual(len(layered),2)
        self.assertEqual({e['plasmid_unit'] for e in layered},{'PU_0001','PU_0002'})
        self.assertTrue(all(e['cluster_relation']=='same_cluster' for e in layered))

    def test_run_multilayer_network_flag_is_opt_in_with_correct_relations(self):
        meta=self.p/'ml_meta.tsv'
        write_tsv(meta,['isolate_id','chromosomal_cluster'],[
            dict(isolate_id='i1',chromosomal_cluster='CC1'),dict(isolate_id='i2',chromosomal_cluster='CC1'),
            dict(isolate_id='i3',chromosomal_cluster='CC2'),dict(isolate_id='i4',chromosomal_cluster='')])
        seq_a='ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC'
        seq_b='TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG'
        rows=[('i1','pA1',seq_a),('i2','pA2',seq_a),('i3','pA3',seq_a),('i4','pA4',seq_a),
              ('i1','pB1',seq_b),('i2','pB2',seq_b)]
        manifest_rows=[]
        for iso,pid,seq in rows:
            path=self.p/f'{pid}.fa';path.write_text(f'>{pid}\n{seq}\n')
            manifest_rows.append(dict(isolate_id=iso,plasmid_id=pid,fasta_path=str(path)))
        manifest=self.p/'ml_manifest.tsv'
        write_tsv(manifest,['isolate_id','plasmid_id','fasta_path'],manifest_rows)
        base=dict(manifest=manifest,metadata=meta,features=None,mode='precomputed',engine='kmer',
            threshold=.05,k=3,min_length=3,sketch_size=100,linkage='complete',annotation_config=None,external_typing=None)
        out_off=self.p/'ml_off'
        run(argparse.Namespace(out=out_off,**base))
        self.assertFalse((out_off/'network.multilayer.graphml').exists())
        self.assertFalse((out_off/'network.multilayer_edges.tsv').exists())
        out_on=self.p/'ml_on'
        run(argparse.Namespace(out=out_on,multilayer_network=True,**base))
        self.assertTrue((out_on/'network.multilayer.graphml').exists())
        layered=read_tsv(out_on/'network.multilayer_edges.tsv')
        pair_i1_i2=[e for e in layered if {e['source'],e['target']}=={'i1','i2'}]
        self.assertEqual(len(pair_i1_i2),2)
        self.assertTrue(all(e['cluster_relation']=='same_cluster' for e in pair_i1_i2))
        pair_i1_i3=[e for e in layered if {e['source'],e['target']}=={'i1','i3'}]
        self.assertEqual(len(pair_i1_i3),1)
        self.assertEqual(pair_i1_i3[0]['cluster_relation'],'cross_cluster')
        pair_i1_i4=[e for e in layered if {e['source'],e['target']}=={'i1','i4'}]
        self.assertEqual(pair_i1_i4[0]['cluster_relation'],'unknown_cluster')

    def test_run_rejects_invalid_containment_bounds_before_reading_inputs(self):
        args=argparse.Namespace(manifest=self.p/'does_not_exist.tsv',metadata=self.p/'also_missing.tsv',
            features=None,mode='precomputed',engine='kmer',out=self.p/'unused_out',threshold=.05,k=3,
            min_length=3,sketch_size=100,linkage='complete',annotation_config=None,external_typing=None,
            containment_min_ratio=0.9,containment_max_ratio=0.5)
        with self.assertRaisesRegex(ValueError,'containment-min-ratio'):
            run(args)
        self.assertFalse((self.p/'unused_out').exists())

    def test_run_writes_margin_and_containment_and_typing_does_not_leak_into_quality(self):
        out1=self.run_fixture('run_no_typing')
        edge_rows=read_tsv(out1/'network.edge_evidence.tsv')
        self.assertEqual(len(edge_rows),1)
        edge=edge_rows[0]
        self.assertAlmostEqual(float(edge['threshold_margin']),.05-float(edge['minimum_distance']))
        self.assertTrue((out1/'containment_candidates.tsv').exists())
        self.assertFalse((out1/'typing_crossreference.tsv').exists())
        self.assertFalse((out1/'network.multilayer.graphml').exists())
        self.assertFalse((out1/'network.multilayer_edges.tsv').exists())
        digest=hashlib.sha256('ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC'.encode()).hexdigest()
        typing_path=self.p/'typing.tsv'
        write_tsv(typing_path,['plasmid_id','sequence_sha256','evidence_source','external_tool','mob_primary_cluster_id'],
                  [dict(plasmid_id='p1',sequence_sha256=digest,evidence_source='lab',external_tool='mob_suite',mob_primary_cluster_id='AA1'),
                   dict(plasmid_id='p2',sequence_sha256=digest,evidence_source='lab',external_tool='mob_suite',mob_primary_cluster_id='AA1')])
        out2=self.run_fixture('run_with_typing',external_typing=typing_path)
        crossref=read_tsv(out2/'typing_crossreference.tsv')
        self.assertEqual({r['plasmid_id'] for r in crossref},{'p1','p2'})
        self.assertEqual(crossref[0]['mob_primary_cluster_id'],'AA1')
        self.assertEqual((out1/'biological_quality.tsv').read_text(),(out2/'biological_quality.tsv').read_text())

    def population_summary_fixture(self):
        meta=self.p/'meta.tsv'
        write_tsv(meta,['isolate_id','organism','location','date'],[
            dict(isolate_id='i1',organism='E. coli',location='SiteA',date='2026-01-01'),
            dict(isolate_id='i2',organism='E. coli',location='SiteA',date='2026-01-05'),
            dict(isolate_id='i3',organism='Klebsiella pneumoniae',location='SiteB',date='2026-02-01')])
        seq_a='ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC'
        seq_b='TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG'
        p1=self.p/'p1.fa';p1.write_text(f'>p1\n{seq_a}\n')
        p2=self.p/'p2.fa';p2.write_text(f'>p2\n{seq_a}\n')
        p3=self.p/'p3.fa';p3.write_text(f'>p3\n{seq_b}\n')
        manifest=self.p/'pop_manifest.tsv'
        write_tsv(manifest,['isolate_id','plasmid_id','fasta_path'],[
            dict(isolate_id='i1',plasmid_id='p1',fasta_path=str(p1)),
            dict(isolate_id='i2',plasmid_id='p2',fasta_path=str(p2)),
            dict(isolate_id='i3',plasmid_id='p3',fasta_path=str(p3))])
        features=self.p/'pop_features.tsv'
        row=dict(plasmid_id='p1',feature_id='f1',start=1,end=10,strand='+',amr_gene='blaTEST',
                 drug_class='beta-lactam',hit_class='strict',annotation_confidence='high',
                 annotation_engine='amrfinder',database_name='amrfinderdb',database_version='2026')
        write_tsv(features,list(row),[row])
        out=self.p/'pop_results'
        args=argparse.Namespace(manifest=manifest,metadata=meta,features=features,mode='precomputed',engine='kmer',
            out=out,threshold=.05,k=3,min_length=3,sketch_size=100,linkage='complete',annotation_config=None)
        run(args)
        return out

    def test_population_summary_rejects_incomplete_or_stale_run(self):
        out=self.population_summary_fixture()
        provenance_path=out/'run_provenance.json'
        record=json.loads(provenance_path.read_text())
        record['status']='running'
        provenance_path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,'requires a completed run'):
            summarize(argparse.Namespace(results=out,out_prefix=None))
        record['status']='complete'
        record['output_sha256']['plasmid_index.tsv']='0'*64
        provenance_path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,'changed or missing'):
            summarize(argparse.Namespace(results=out,out_prefix=None))

    def test_population_summary_aggregates_expected_counts(self):
        out=self.population_summary_fixture()
        summarize(argparse.Namespace(results=out,out_prefix=None))
        pu_rows=read_tsv(out/'population_summary.pu_level.tsv')
        self.assertEqual(len(pu_rows),2)
        pu_rows.sort(key=lambda r:-int(r['n_plasmids']))
        big,small=pu_rows
        self.assertEqual(big['n_plasmids'],'2')
        self.assertEqual(set(big['isolates'].split(';')),{'i1','i2'})
        self.assertEqual(big['resistance_genes'],'blaTEST')
        self.assertEqual(big['n_resistance_genes'],'1')
        self.assertEqual(small['n_plasmids'],'1')
        self.assertEqual(small['isolates'],'i3')
        self.assertEqual(small['resistance_genes'],'')
        dims=read_tsv(out/'population_summary.metadata_dimension.tsv')
        organism_rows={r['value']:r for r in dims if r['dimension']=='organism'}
        self.assertEqual(organism_rows['E. coli']['n_isolates'],'2')
        self.assertEqual(organism_rows['Klebsiella pneumoniae']['n_isolates'],'1')
        self.assertNotIn('chromosomal_cluster',{r['dimension'] for r in dims})


if __name__=='__main__':unittest.main()
