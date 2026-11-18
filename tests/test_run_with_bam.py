import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
SPEC = importlib.util.spec_from_file_location('run_with_bam', ROOT / 'tools/run_with_bam.py')
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)
from reference_path import populate_reference_path, read_reference_sequences


class TemporaryBamTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.samtools = shutil.which('samtools')
        cls.tool = Path(os.environ.get('TOOL', ROOT / 'build/isamtools')).resolve()
        if not cls.samtools or not cls.tool.is_file():
            raise RuntimeError('samtools and built isamtools are required')
        cls.fixture = tempfile.TemporaryDirectory(prefix='isam-with-bam-test-')
        cls.bundle = Path(cls.fixture.name)
        cls.reference = cls.bundle / 'reference.fa'
        cls.reference.write_text('>chr1\n' + 'ACGT' * 250 + '\n')
        sam = cls.bundle / 'source.sam'
        sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
                       '@RG\tID:movie\tSM:synthetic\n' + ''.join(
                           '\t'.join([f'read{i}', str(flag), 'chr1' if flag != 4 else '*',
                                      str(i * 100 + 1) if flag != 4 else '0', '60',
                                      '8M' if flag != 4 else '*', '*', '0', '0',
                                      'ACGTACGT', 'IIIIIIII', 'RG:Z:movie', 'HP:i:1',
                                      'PS:i:17', 'MM:Z:C+m,0;', 'ML:B:C,200']) + '\n'
                           for i, flag in enumerate((0, 0, 2048, 4))))
        commands = [[cls.samtools, 'faidx', str(cls.reference)],
                    [cls.samtools, 'view', '-b', '-o', str(cls.bundle / 'source.bam'), str(sam)],
                    [cls.samtools, 'view', '-C', '-T', str(cls.reference), '-o', str(cls.bundle / 'canon.cram'), str(sam)],
                    [cls.samtools, 'index', str(cls.bundle / 'canon.cram')],
                    [str(cls.tool), 'index', '-T', str(cls.reference), str(cls.bundle / 'canon.cram')],
                    [str(cls.tool), 'import', '--fast-verify', '--store', str(cls.bundle / 'canon.cram'),
                     '--store-index', str(cls.bundle / 'canon.ssam'), '--aln', str(cls.bundle / 'hifi.hg38.isam'),
                     str(cls.bundle / 'source.bam')],
                    [str(cls.tool), 'sort', '--all', str(cls.bundle / 'hifi.hg38.isam')]]
        with (cls.bundle / 'setup.log').open('wb') as log:
            for command in commands:
                subprocess.run(command, check=True, stdout=log, stderr=log)

    @classmethod
    def tearDownClass(cls):
        cls.fixture.cleanup()

    def arguments(self, base, **overrides):
        values = dict(bundle=self.bundle, view='hifi.hg38', reference=self.reference, reference_path=None,
                      scratch=base, isamtools=str(self.tool), samtools=self.samtools,
                      threads=3, decode_threads=0, cache_mb=16, stage='none', stage_reference=False, strategy='coordinate',
                      decode_lookahead=None, cache_selected=False, write_index=False,
                      region=None, expected_records=None, report=base / 'report.json',
                      command=[self.samtools, 'view', '-h', '-o', str(base / 'consumer.sam'), '{bam}'])
        values.update(overrides)
        return argparse.Namespace(**values)

    def verify_result(self, base):
        report = json.loads((base / 'report.json').read_text())
        self.assertEqual(report['records'], 4)
        self.assertTrue(report['local_bulk_data_removed'])
        self.assertFalse(Path(report['scratch']).exists())
        records = [line for line in (base / 'consumer.sam').read_text().splitlines() if not line.startswith('@')]
        self.assertEqual(len(records), 4)
        for record in records:
            fields = record.split('\t')
            self.assertEqual(fields[9:11], ['ACGTACGT', 'IIIIIIII'])
            for tag in ['RG:Z:movie', 'HP:i:1', 'PS:i:17', 'MM:Z:C+m,0;', 'ML:B:C,200']:
                self.assertIn(tag, fields[11:])
        self.assertEqual(sum(bool(int(record.split('\t')[1]) & 4) for record in records), 1)

    def test_independent_coordinate_decoders_preserve_tags_and_cleanup(self):
        with tempfile.TemporaryDirectory(prefix='isam-with-decoder-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, decode_threads=4)), 0)
            self.verify_result(base)
            self.assertEqual(json.loads((base / 'report.json').read_text())['decode_threads'], 4)

    def test_decode_workers_require_coordinate_isam_before_staging(self):
        for overrides in ({'strategy': 'store-sort'}, {'view': 'canonical', 'expected_records': 4}):
            with tempfile.TemporaryDirectory(prefix='isam-invalid-decoders-') as directory:
                base = Path(directory)
                with self.assertRaisesRegex(RuntimeError, 'requires --strategy coordinate'):
                    helper.execute(self.arguments(base, decode_threads=4, **overrides))
                self.assertFalse((base / 'report.json').exists())

    def test_selected_cache_and_inline_index_reach_unmodified_consumer(self):
        with tempfile.TemporaryDirectory(prefix='isam-selected-consumer-') as directory:
            base = Path(directory)
            commands = []
            original = helper.checked
            def checked(command, env, capture=False):
                commands.append(command)
                return original(command, env, capture)
            with patch.object(Path, 'is_relative_to', return_value=True), patch.object(helper, 'checked', checked):
                self.assertEqual(helper.execute(self.arguments(
                    base, decode_threads=4, decode_lookahead=17,
                    cache_selected=True, write_index=True)), 0)
            self.verify_result(base)
            report = json.loads((base / 'report.json').read_text())
            self.assertTrue(report['index_included_in_conversion'])
            self.assertEqual(report['index_seconds'], 0)
            self.assertFalse(any(c[:2] == [self.samtools, 'index'] for c in commands))

    def test_new_cache_options_reject_incompatible_routes_before_staging(self):
        for overrides in ({'cache_selected': True}, {'decode_lookahead': 17},
                          {'write_index': True, 'strategy': 'store-sort'},
                          {'write_index': True, 'view': 'canonical', 'expected_records': 4}):
            with tempfile.TemporaryDirectory(prefix='isam-invalid-cache-') as directory:
                base = Path(directory)
                with self.assertRaisesRegex(RuntimeError, 'require'):
                    helper.execute(self.arguments(base, **overrides))
                self.assertFalse((base / 'report.json').exists())

    def test_m5_coordinate_decoders_do_not_need_the_original_fasta(self):
        with tempfile.TemporaryDirectory(prefix='isam-with-m5-decoders-') as directory:
            base = Path(directory)
            reference_cache = base / 'm5'
            populate_reference_path(self.reference,
                                    read_reference_sequences(self.bundle / 'canon.cram', self.samtools),
                                    reference_cache)
            hidden = self.reference.with_suffix('.hidden')
            self.reference.rename(hidden)
            try:
                with patch.object(Path, 'is_relative_to', return_value=True):
                    self.assertEqual(helper.execute(self.arguments(
                        base, reference=None, reference_path=reference_cache,
                        decode_threads=4, stage_reference=True)), 0)
                self.verify_result(base)
            finally:
                hidden.rename(self.reference)

    def test_all_input_and_reference_staging_preserves_tags_and_cleans(self):
        before = {name: (self.bundle / name).read_bytes()
                  for name in ['canon.cram', 'canon.ssam', 'hifi.hg38.isam']}
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, stage='all', stage_reference=True)), 0)
            self.verify_result(base)
        for name, contents in before.items():
            self.assertEqual((self.bundle / name).read_bytes(), contents)

    def test_indexes_only_staging(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True), \
                    patch.object(helper, 'checked', wraps=helper.checked) as checked:
                self.assertEqual(helper.execute(self.arguments(base, stage='indexes')), 0)
            self.verify_result(base)
            report = json.loads((base / 'report.json').read_text())
            count_command = next(call.args[0] for call in checked.call_args_list
                                 if '-c' in call.args[0])
            self.assertEqual(count_command[-1], str(Path(report['scratch']) / 'inputs/hifi.hg38.isam'))
            self.assertEqual(count_command[count_command.index('--store-index') + 1],
                             str(Path(report['scratch']) / 'inputs/canon.ssam'))

    def test_canonical_keeps_unmapped_and_all_alignment_records(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, view='canonical', expected_records=4)), 0)
            self.verify_result(base)

    def test_store_sort_keeps_all_payloads_and_cleans_sort_intermediates(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, strategy='store-sort', stage='indexes')), 0)
            self.verify_result(base)
            self.assertEqual(json.loads((base / 'report.json').read_text())['strategy'], 'store-sort')

    def test_store_sort_producer_failure_prevents_consumer_and_cleans(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            # Counts are supplied, so this fails only in the producer stage.
            with patch.object(Path, 'is_relative_to', return_value=True):
                with self.assertRaisesRegex(RuntimeError, 'store-sort failed: producer=1'):
                    helper.execute(self.arguments(base, strategy='store-sort', isamtools='/bin/false', expected_records=4))
            self.assertFalse((base / 'consumer.sam').exists())
            report = json.loads((base / 'report.json').read_text())
            self.assertTrue(report['local_bulk_data_removed'])
            self.assertFalse(Path(report['scratch']).exists())

    def test_store_sort_cli_refuses_region(self):
        result = subprocess.run([sys.executable, str(ROOT / 'tools/run_with_bam.py'),
                                 '--bundle', str(self.bundle), '--view', 'hifi.hg38',
                                 '--reference', str(self.reference), '--scratch', self.fixture.name,
                                 '--strategy', 'store-sort', '--region', 'chr1:1-10',
                                 '--', 'consumer', '{bam}'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('store-sort requires a whole iSAM view', result.stderr)

    def test_bad_count_prevents_consumer_and_cleans(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                with self.assertRaisesRegex(RuntimeError, 'record count mismatch'):
                    helper.execute(self.arguments(base, expected_records=5))
            self.assertFalse((base / 'consumer.sam').exists())
            report = json.loads((base / 'report.json').read_text())
            self.assertTrue(report['local_bulk_data_removed'])
            self.assertFalse(Path(report['scratch']).exists())

    def test_failed_consumer_returns_its_code_and_cleans(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            command = [sys.executable, '-c', 'import sys; sys.exit(7)', '{bam}']
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, command=command)), 7)
            report = json.loads((base / 'report.json').read_text())
            self.assertTrue(report['local_bulk_data_removed'])
            self.assertFalse(Path(report['scratch']).exists())

    def test_shared_storage_is_refused(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            with self.assertRaisesRegex(RuntimeError, '/scratch-node'):
                helper.execute(self.arguments(Path(directory)))

    def test_nonstaged_shared_indexes_are_refused_and_cleaned(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', autospec=True,
                              side_effect=lambda path, other: path == base):
                with self.assertRaisesRegex(RuntimeError, 'use --stage indexes'):
                    helper.execute(self.arguments(base, stage='none'))
            report = json.loads((base / 'report.json').read_text())
            self.assertTrue(report['local_bulk_data_removed'])
            self.assertFalse(Path(report['scratch']).exists())

    def test_m5_staging_preserves_tags_sources_and_cleanup_for_both_strategies(self):
        with tempfile.TemporaryDirectory(prefix='isam-m5-helper-test-') as directory:
            base = Path(directory)
            cache = base / 'cache'
            populate_reference_path(self.reference,
                                    read_reference_sequences(self.bundle / 'canon.cram', self.samtools), cache)
            before = {path.name: path.read_bytes() for path in cache.iterdir()}
            for strategy in ('coordinate', 'store-sort'):
                output = base / strategy
                output.mkdir()
                with patch.object(Path, 'is_relative_to', return_value=True):
                    self.assertEqual(helper.execute(self.arguments(
                        output, reference=None, reference_path=cache, strategy=strategy,
                        stage='indexes', stage_reference=True)), 0)
                self.verify_result(output)
                report = json.loads((output / 'report.json').read_text())
                self.assertEqual(report['reference_mode'], 'm5')
                self.assertGreater(report['conversion_children_cpu_seconds'], 0)
                self.assertEqual(report['conversion_children_cpu_seconds'],
                                 report['conversion_children_user_seconds'] +
                                 report['conversion_children_system_seconds'])
            self.assertEqual({path.name: path.read_bytes() for path in cache.iterdir()}, before)

    def test_m5_canonical_runs_without_the_original_fasta(self):
        with tempfile.TemporaryDirectory(prefix='isam-m5-canonical-test-') as directory:
            base = Path(directory)
            cache = base / 'cache'
            populate_reference_path(self.reference,
                                    read_reference_sequences(self.bundle / 'canon.cram', self.samtools), cache)
            hidden = self.reference.with_name('hidden-reference.fa')
            fai = Path(str(self.reference) + '.fai')
            hidden_fai = Path(str(hidden) + '.fai')
            self.reference.rename(hidden)
            fai.rename(hidden_fai)
            try:
                with patch.object(Path, 'is_relative_to', return_value=True):
                    self.assertEqual(helper.execute(self.arguments(
                        base, reference=None, reference_path=cache, view='canonical', expected_records=4)), 0)
                self.verify_result(base)
            finally:
                hidden.rename(self.reference)
                hidden_fai.rename(fai)

    def test_missing_m5_prevents_consumer_and_preserves_the_cache(self):
        with tempfile.TemporaryDirectory(prefix='isam-m5-missing-test-') as directory:
            base = Path(directory)
            cache = base / 'cache'
            cache.mkdir()
            with patch.object(Path, 'is_relative_to', return_value=True):
                with self.assertRaisesRegex(RuntimeError, 'missing/wrong-length M5 reference'):
                    helper.execute(self.arguments(base, reference=None, reference_path=cache, stage='indexes'))
            report = json.loads((base / 'report.json').read_text())
            self.assertTrue(report['local_bulk_data_removed'])
            self.assertFalse((base / 'consumer.sam').exists())
            self.assertTrue(cache.is_dir())
            self.assertEqual(list(cache.iterdir()), [])

    def test_reference_builder_refuses_existing_data_and_removes_own_failed_output(self):
        with tempfile.TemporaryDirectory(prefix='isam-m5-preparation-test-') as directory:
            base = Path(directory)
            existing = base / 'existing'
            existing.mkdir()
            (existing / 'user-file').write_bytes(b'keep')
            sequences = read_reference_sequences(self.bundle / 'canon.cram', self.samtools)
            with self.assertRaises(FileExistsError):
                populate_reference_path(self.reference, sequences, existing)
            self.assertEqual((existing / 'user-file').read_bytes(), b'keep')
            wrong = [dict(sequences[0], M5='0' * 32)]
            with self.assertRaisesRegex(RuntimeError, 'canonical reference M5'):
                populate_reference_path(self.reference, wrong, base / 'failed')
            self.assertFalse((base / 'failed').exists())

    def test_reference_builder_cli_prepares_verified_files_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory(prefix='isam-m5-cli-test-') as directory:
            cache = Path(directory) / 'cache'
            command = [sys.executable, str(ROOT / 'tools/reference_path.py'),
                       '--cram', str(self.bundle / 'canon.cram'),
                       '--reference', str(self.reference), '--output', str(cache),
                       '--samtools', self.samtools]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(result.stdout)
            self.assertEqual(report['contigs_verified'], 1)
            self.assertEqual(report['bytes'], 1000)
            before = {path.name: path.read_bytes() for path in cache.iterdir()}
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('reference preparation failed', result.stderr)
            self.assertEqual({path.name: path.read_bytes() for path in cache.iterdir()}, before)

    def test_reference_path_special_characters_fail_before_consumer_or_population(self):
        with tempfile.TemporaryDirectory(prefix='isam-m5-path-test-') as directory:
            base = Path(directory)
            sequences = read_reference_sequences(self.bundle / 'canon.cram', self.samtools)
            for name in ('cache:second', 'cache%s'):
                with self.subTest(name=name):
                    cache = base / name
                    with self.assertRaisesRegex(RuntimeError, 'cannot contain : or %'):
                        populate_reference_path(self.reference, sequences, cache)
                    self.assertFalse(cache.exists())
                    cache.mkdir()
                    output = base / ('output-' + name)
                    output.mkdir()
                    with patch.object(Path, 'is_relative_to', return_value=True):
                        with self.assertRaisesRegex(RuntimeError, 'without : or %'):
                            helper.execute(self.arguments(output, reference=None,
                                                          reference_path=cache, stage='indexes'))
                    self.assertFalse((output / 'consumer.sam').exists())
                    report = json.loads((output / 'report.json').read_text())
                    self.assertTrue(report['local_bulk_data_removed'])
                    self.assertFalse(Path(report['scratch']).exists())
                    self.assertTrue(cache.is_dir())

    def test_reference_cli_refuses_both_reference_forms(self):
        result = subprocess.run([sys.executable, str(ROOT / 'tools/run_with_bam.py'),
                                 '--bundle', str(self.bundle), '--view', 'hifi.hg38',
                                 '--reference', str(self.reference), '--reference-path', str(self.bundle),
                                 '--scratch', self.fixture.name, '--', 'consumer', '{bam}'],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not allowed with argument', result.stderr)


class LocalIndexGuardTests(unittest.TestCase):
    def setUp(self):
        self.fixture = tempfile.TemporaryDirectory(prefix='isam-helper-path-test-')
        self.addCleanup(self.fixture.cleanup)
        self.base = Path(self.fixture.name)
        self.local = self.base / 'local'
        self.shared = self.base / 'shared'
        self.work = self.base / 'work'
        for path in (self.local, self.shared, self.work):
            path.mkdir()
        for name in ('canon.cram', 'reference.fa', 'reference.fa.fai'):
            (self.shared / name).write_bytes(b'tiny-' + name.encode())
        (self.local / 'canon.cram').symlink_to(self.shared / 'canon.cram')
        for name in ('canon.ssam', 'hifi.hg38.isam'):
            (self.local / name).write_bytes(b'tiny-local-index')
            (self.shared / name).write_bytes(b'tiny-shared-index')
        self.args = argparse.Namespace(
            bundle=self.local, view='hifi.hg38', reference=self.shared / 'reference.fa', reference_path=None,
            stage='none', region=None, expected_records=None, stage_reference=False,
            strategy='coordinate', threads=1, decode_threads=0, decode_lookahead=None,
            cache_selected=False, write_index=False, cache_mb=1, isamtools='unused', samtools='unused')

        original = Path.is_relative_to
        self.locality = patch.object(
            Path, 'is_relative_to', autospec=True,
            side_effect=lambda path, root: original(path, self.local) if root == Path('/scratch-node')
            else original(path, root))

    def reject_shared_symlink(self, name):
        path = self.local / name
        path.unlink()
        path.symlink_to(self.shared / name)
        with self.locality, patch.object(helper, 'checked', side_effect=RuntimeError('count reached')) as checked:
            with self.assertRaisesRegex(RuntimeError, 'use --stage indexes'):
                helper.materialize(self.args, self.work, {'staging': []})
            checked.assert_not_called()
        self.assertFalse((self.work / 'inputs').exists())
        self.assertTrue(path.is_symlink())
        self.assertTrue((self.shared / name).is_file())

    def test_nonstaged_ssam_symlink_to_shared_storage_is_rejected_before_count(self):
        self.reject_shared_symlink('canon.ssam')

    def test_nonstaged_isam_symlink_to_shared_storage_is_rejected_before_count(self):
        self.reject_shared_symlink('hifi.hg38.isam')

    def test_local_lmdb_allows_shared_cram_and_reference(self):
        with self.locality, patch.object(helper, 'checked', side_effect=RuntimeError('count reached')) as checked:
            with self.assertRaisesRegex(RuntimeError, 'count reached'):
                helper.materialize(self.args, self.work, {'staging': []})
        command = checked.call_args.args[0]
        self.assertIn('-c', command)
        self.assertEqual(command[command.index('--store-index') + 1], str(self.local / 'canon.ssam'))
        self.assertEqual(command[-1], str(self.local / 'hifi.hg38.isam'))
        self.assertTrue((self.local / 'canon.cram').is_symlink())

    def test_auto_strategy_matches_scope_and_preserves_explicit_choice(self):
        for view, region, requested, expected in (
                ('hifi.hg38', None, 'auto', 'store-sort'),
                ('nonhifi.chm13', None, 'auto', 'store-sort'),
                ('hifi.hg38', 'chr1:1-200', 'auto', 'coordinate'),
                ('canonical', None, 'auto', 'coordinate'),
                ('canonical', 'chr1:1-200', 'auto', 'coordinate'),
                ('hifi.hg38', None, 'coordinate', 'coordinate'),
                ('hifi.hg38', None, 'store-sort', 'store-sort')):
            with self.subTest(view=view, region=region, requested=requested):
                self.assertEqual(helper.resolve_strategy(view, region, requested), expected)


class AutomaticStrategyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        TemporaryBamTests.setUpClass.__func__(cls)

    @classmethod
    def tearDownClass(cls):
        TemporaryBamTests.tearDownClass.__func__(cls)

    arguments = TemporaryBamTests.arguments
    verify_result = TemporaryBamTests.verify_result

    def test_auto_whole_view_materializes_sorts_indexes_and_cleans(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-auto-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, strategy='auto', stage='indexes')), 0)
            self.verify_result(base)
            report = json.loads((base / 'report.json').read_text())
            self.assertEqual(report['requested_strategy'], 'auto')
            self.assertEqual(report['strategy'], 'store-sort')

    def test_auto_region_uses_indexed_coordinate_output(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-auto-region-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True), \
                    patch.object(helper, 'checked', wraps=helper.checked) as checked:
                self.assertEqual(helper.execute(self.arguments(base, strategy='auto', stage='indexes',
                                                               region='chr1:1-200')), 0)
            report = json.loads((base / 'report.json').read_text())
            self.assertEqual(report['records'], 2)
            self.assertEqual(report['strategy'], 'coordinate')
            self.assertTrue(report['local_bulk_data_removed'])
            self.assertFalse(Path(report['scratch']).exists())
            producer = next(call.args[0] for call in checked.call_args_list
                            if '--order' in call.args[0])
            self.assertEqual(producer[producer.index('--order') + 1], 'coordinate')
            self.assertEqual(producer[-1], 'chr1:1-200')

    def test_auto_canonical_retains_all_alignment_records(self):
        with tempfile.TemporaryDirectory(prefix='isam-consumer-auto-canon-test-') as directory:
            base = Path(directory)
            with patch.object(Path, 'is_relative_to', return_value=True):
                self.assertEqual(helper.execute(self.arguments(base, view='canonical', strategy='auto',
                                                               expected_records=4)), 0)
            self.verify_result(base)
            self.assertEqual(json.loads((base / 'report.json').read_text())['strategy'], 'coordinate')

    def test_cli_default_requests_automatic_selection(self):
        command = ['run_with_bam.py', '--bundle', str(self.bundle), '--view', 'hifi.hg38',
                   '--reference', str(self.reference), '--scratch', self.fixture.name,
                   '--', 'consumer', '{bam}']
        with patch.object(sys, 'argv', command), patch.object(helper, 'execute', return_value=0) as execute:
            self.assertEqual(helper.main(), 0)
        self.assertEqual(execute.call_args.args[0].strategy, 'auto')


if __name__ == '__main__':
    unittest.main()
