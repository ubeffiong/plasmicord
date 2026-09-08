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
from python.contracts import write_tsv
from python.calibration import calibrate
from python.cli import aggregate


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


if __name__=='__main__':unittest.main()
