"""Coordinate BAMs stay ordered under cache eviction and parallel decoding."""
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class CoordinateDecodeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tool = str(Path(os.environ.get('TOOL', ROOT / 'build/isamtools')).resolve())
        cls.samtools = shutil.which('samtools')
        if not cls.samtools or not Path(cls.tool).is_file():
            raise RuntimeError('samtools and built isamtools are required')
        cls.fixture = tempfile.TemporaryDirectory(prefix='isam-coordinate-decode-')
        cls.base = Path(cls.fixture.name)
        cls.reference = cls.base / 'reference.fa'
        cls.reference.write_text('>chr1\n' + 'ACGT' * 3000 + '\n>chr2\n' + 'ACGT' * 3000 + '\n')
        header = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:12000\n@SQ\tSN:chr2\tLN:12000\n' \
                 '@RG\tID:canon\tSM:test\n@RG\tID:view\tSM:test\n'
        view_header = header.replace('@SQ\tSN:chr1\tLN:12000\n@SQ\tSN:chr2\tLN:12000',
                                     '@SQ\tSN:chr2\tLN:12000\n@SQ\tSN:chr1\tLN:12000')
        sequence, quality = 'ACGTTGCA' * 250, 'IHGFEDCB' * 250

        def row(name, flag, contig, position, rg, hp):
            seq = 'TGCAACGT' * 250 if flag & 16 else sequence
            qual = quality[::-1] if flag & 16 else quality
            return '\t'.join([name, str(flag), contig, str(position), '60',
                              '2000M' if contig != '*' else '*', '*', '0', '0', seq, qual,
                              'RG:Z:' + rg, 'HP:i:' + str(hp), 'PS:i:17', 'MM:Z:C+m,0;',
                              'ML:B:C,200', 'zm:i:101', 'rq:f:0.999']) + '\n'

        canonical, view = [header], [view_header]
        for i in range(600):
            name = 'read' + str(i)
            canonical.append(row(name, 0, 'chr1', 1 + 10 * i, 'canon', 1))
            view.append(row(name, 16 if i % 2 else 0, 'chr1', 6001 - 10 * i, 'view', 2))
            if i % 13 == 0:
                view.append(row(name, 2048, 'chr2', 1 + 10 * i, 'view', 2))
            if i % 17 == 0:
                view.append(row(name, 256, 'chr1', 6001 - 10 * i, 'view', 2))
        for i in range(3):
            canonical.append(row('unmapped' + str(i), 4, '*', 0, 'canon', 1))
            view.append(row('unmapped' + str(i), 4, '*', 0, 'view', 2))
        (cls.base / 'canonical.sam').write_text(''.join(canonical))
        (cls.base / 'view.sam').write_text(''.join(view))
        cls.checked([cls.samtools, 'faidx', str(cls.reference)])
        cls.store, cls.index, cls.view = [cls.base / name for name in ('canonical.cram', 'canonical.ssam', 'view.isam')]
        cls.checked([cls.samtools, 'view', '-C', '-T', str(cls.reference),
                     '--output-fmt-option', 'seqs_per_slice=12',
                     '--output-fmt-option', 'slices_per_container=1',
                     '-o', str(cls.store), str(cls.base / 'canonical.sam')])
        cls.checked([cls.samtools, 'index', str(cls.store)])
        cls.checked([cls.samtools, 'view', '-b', '-o', str(cls.base / 'view.bam'), str(cls.base / 'view.sam')])
        cls.checked([cls.tool, 'index', '-T', str(cls.reference), '-o', str(cls.index), str(cls.store)])
        cls.checked([cls.tool, 'import', '--fast-verify', '--store', str(cls.store),
                     '--store-index', str(cls.index), '--aln', str(cls.view), str(cls.base / 'view.bam')])
        cls.checked([cls.tool, 'sort', '--all', str(cls.view)])

    @classmethod
    def tearDownClass(cls):
        cls.fixture.cleanup()

    @staticmethod
    def checked(command, env=None):
        result = subprocess.run(command, capture_output=True, text=True, env=env, timeout=60)
        if result.returncode:
            raise RuntimeError(f'{command}: {result.stderr}')
        return result

    def emit(self, workers=None, options=(), cache_mb=1):
        command = [self.tool, 'view', '-h', '-T', str(self.reference), '--order', 'coordinate',
                   '--store', str(self.store), '--store-index', str(self.index)]
        if workers is not None:
            command += ['--decode-threads', str(workers)]
        env = dict(os.environ, ISAM_SLICE_CACHE_MB=str(cache_mb), ISAM_FETCH_STATS='1')
        return self.checked(command + list(options) + [str(self.view)], env)

    def test_parallel_metadata_payloads_and_order(self):
        expected = self.emit().stdout
        for workers in (1, 2, 4):
            for cache_mb in (1, 8):
                with self.subTest(workers=workers, cache_mb=cache_mb):
                    result = self.emit(workers, cache_mb=cache_mb)
                    self.assertEqual(result.stdout, expected)
                    self.assertGreater(int(re.search(r'miss/decode=(\d+)', result.stderr).group(1)), 1)

    @staticmethod
    def pipeline_stats(result):
        line = re.search(r'pipeline-stats: (.*)', result.stderr)
        if line is None:
            raise AssertionError('missing asynchronous pipeline statistics')
        return {key: float(value) for key, value in re.findall(r'(\w+)=([0-9.]+)', line.group(1))}

    def test_pipeline_pending_results_are_bounded_and_drained(self):
        expected = self.emit().stdout
        for workers in (1, 2, 4, 64):
            with self.subTest(workers=workers):
                result = self.emit(workers)
                self.assertEqual(result.stdout, expected)
                stats = self.pipeline_stats(result)
                self.assertGreater(stats['submitted'], 1)
                self.assertEqual(stats['submitted'], stats['completed'])
                self.assertEqual(stats['pending_limit'], workers * 2)
                self.assertGreater(stats['peak_pending'], 0)
                self.assertLessEqual(stats['peak_pending'], stats['pending_limit'])

    def test_filtered_window_does_not_submit_unused_containers(self):
        result = self.emit(4, ('-q', '100'))
        self.assertEqual(result.stdout, self.emit(options=('-q', '100')).stdout)
        stats = self.pipeline_stats(result)
        self.assertEqual(stats['submitted'], 0)
        self.assertEqual(stats['completed'], 0)
        self.assertEqual(stats['peak_pending'], 0)

    def test_oversized_container_survives_speculative_decodes(self):
        store, index, view = [self.base / name for name in
                              ('oversized.cram', 'oversized.ssam', 'oversized.isam')]
        self.checked([self.samtools, 'view', '-C', '-T', str(self.reference),
                      '--output-fmt-option', 'seqs_per_slice=600',
                      '--output-fmt-option', 'slices_per_container=1',
                      '-o', str(store), str(self.base / 'canonical.sam')])
        self.checked([self.samtools, 'index', str(store)])
        self.checked([self.tool, 'index', '-T', str(self.reference), '-o', str(index), str(store)])
        self.checked([self.tool, 'import', '--fast-verify', '--store', str(store),
                      '--store-index', str(index), '--aln', str(view), str(self.base / 'view.bam')])
        self.checked([self.tool, 'sort', '--all', str(view)])
        expected = self.emit().stdout
        for workers in (1, 2, 4):
            command = [self.tool, 'view', '-h', '-T', str(self.reference),
                       '--decode-threads', str(workers), '--store', str(store),
                       '--store-index', str(index), str(view)]
            result = self.checked(command, dict(os.environ, ISAM_SLICE_CACHE_MB='1', ISAM_FETCH_STATS='1'))
            self.assertEqual(result.stdout, expected)
            peak = int(re.search(r'cache_peak_bytes=(\d+)', result.stderr).group(1))
            self.assertGreater(peak, 1024 * 1024)
            stats = self.pipeline_stats(result)
            self.assertEqual(stats['submitted'], stats['completed'])

    def test_filters_and_single_region(self):
        for options in (('-F', '256'), ('-f', '2048'), ('-r', 'view'), ('-d', 'HP:2'),
                        ('-x', 'MM,ML'), ('-q', '100')):
            with self.subTest(options=options):
                self.assertEqual(self.emit(4, options).stdout, self.emit(options=options).stdout)
        for region in ('chr1:1-3000', 'chr2:1-3000', '*'):
            command = [self.tool, 'view', '-h', '-T', str(self.reference), '--store', str(self.store),
                       '--store-index', str(self.index)]
            serial = self.checked(command + [str(self.view), region]).stdout
            parallel = self.checked(command + ['--decode-threads', '4', str(self.view), region]).stdout
            self.assertEqual(parallel, serial)

    def test_direct_bam_is_indexable_without_sort(self):
        output = self.base / 'direct.bam'
        result = self.emit(4, ('-1', '-o', str(output)))
        self.assertEqual(result.stdout, '')
        self.checked([self.samtools, 'quickcheck', '-v', str(output)])
        self.checked([self.samtools, 'index', str(output)])
        self.assertEqual(int(self.checked([self.samtools, 'view', '-c', str(output)]).stdout), 686)

    @staticmethod
    def selection_stats(result):
        line = re.search(r'selection-stats: (.*)', result.stderr)
        if line is None:
            raise AssertionError('missing selective-cache statistics')
        return {key: float(value) for key, value in re.findall(r'(\w+)=([0-9.]+)', line.group(1))}

    def test_selected_cache_preserves_repeated_payloads_and_releases_final_uses(self):
        expected = self.emit().stdout
        for workers, window in ((1, 1), (2, 17), (4, 512), (4, 4096)):
            with self.subTest(workers=workers, window=window):
                result = self.emit(workers, ('--cache-selected', '--decode-lookahead', str(window)))
                self.assertEqual(result.stdout, expected)
                stats = self.selection_stats(result)
                self.assertEqual(stats['requests'], 686)
                self.assertEqual(stats['remaining'], 0)
                self.assertEqual(stats['released_records'], 603)
                pending = self.pipeline_stats(result)
                self.assertEqual(pending['submitted'], pending['completed'])
                self.assertLessEqual(pending['peak_pending'], workers * 2)

    def test_selected_cache_with_metadata_and_payload_filters(self):
        for options in (('-F', '2304'), ('-f', '2048'), ('-q', '100'),
                        ('-r', 'view'), ('-d', 'HP:1'), ('-x', 'MM,ML'),
                        ('--subsample', '0.3', '--subsample-seed', '17')):
            with self.subTest(options=options):
                result = self.emit(4, ('--cache-selected', *options))
                self.assertEqual(result.stdout, self.emit(options=options).stdout)
                self.assertEqual(self.selection_stats(result)['remaining'], 0)
                if options == ('-q', '100'):
                    self.assertEqual(self.selection_stats(result)['requests'], 0)
                    self.assertEqual(self.pipeline_stats(result)['submitted'], 0)

    def test_selected_cache_region_and_unmapped_payloads(self):
        for region in ('chr1:1-3000', 'chr2:1-3000', '*', 'chr1:11000-12000'):
            command = [self.tool, 'view', '-h', '-T', str(self.reference),
                       '--store', str(self.store), '--store-index', str(self.index)]
            expected = self.checked(command + [str(self.view), region]).stdout
            result = self.checked(command + ['--decode-threads', '4', '--cache-selected',
                                            str(self.view), region],
                                  dict(os.environ, ISAM_SLICE_CACHE_MB='1', ISAM_FETCH_STATS='1'))
            self.assertEqual(result.stdout, expected)
            self.assertEqual(self.selection_stats(result)['remaining'], 0)

    def test_selected_cache_bam_remains_indexable_without_sort(self):
        output = self.base / 'selected.bam'
        self.emit(4, ('--cache-selected', '--decode-lookahead', '17', '-1', '-o', str(output)))
        self.checked([self.samtools, 'quickcheck', '-v', str(output)])
        self.checked([self.samtools, 'index', str(output)])
        self.assertEqual(int(self.checked([self.samtools, 'view', '-c', str(output)]).stdout), 686)

    def test_inline_bai_matches_independent_index_and_payload(self):
        for workers, selected in ((None, False), (4, False), (4, True)):
            with self.subTest(workers=workers, selected=selected):
                output = self.base / 'inline.bam'
                options = ['--write-index', '-1', '-@', '3', '-o', str(output)]
                if selected:
                    options += ['--cache-selected', '--decode-lookahead', '17']
                self.emit(workers, options)
                index = Path(str(output) + '.bai')
                self.assertTrue(index.is_file())
                self.checked([self.samtools, 'quickcheck', '-v', str(output)])
                counts = self.checked([self.samtools, 'idxstats', str(output)]).stdout
                self.assertEqual(sum(int(row.split('\t')[2]) + int(row.split('\t')[3])
                                     for row in counts.splitlines()), 686)
                probe = self.checked([self.samtools, 'view', str(output), 'chr1:1-3000']).stdout
                self.checked([self.samtools, 'index', str(output)])
                self.assertEqual(self.checked([self.samtools, 'idxstats', str(output)]).stdout, counts)
                self.assertEqual(self.checked([self.samtools, 'view', str(output), 'chr1:1-3000']).stdout, probe)

    def test_inline_index_empty_filtered_output(self):
        output = self.base / 'inline-empty.bam'
        self.emit(4, ('--cache-selected', '--write-index', '-q', '100', '-1', '-o', str(output)))
        self.assertTrue(Path(str(output) + '.bai').is_file())
        self.assertEqual(self.checked([self.samtools, 'view', '-c', str(output), 'chr1']).stdout.strip(), '0')

    def test_invalid_combinations_fail(self):
        for options in (('--decode-threads', '0'), ('--decode-threads', '65'),
                        ('--decode-threads', '4294967297'),
                        ('--decode-threads', '4', '--order', 'store'),
                        ('--decode-threads', '4', '-c'), ('--decode-threads', '4', '-p'),
                        ('--cache-selected',), ('--decode-lookahead', '4096'),
                        ('--decode-threads', '4', '--decode-lookahead', '0'),
                        ('--decode-threads', '4', '--decode-lookahead', '65537'),
                        ('--write-index',), ('--write-index', '-b'),
                        ('--write-index', '-b', '-o', '/dev/null'),
                        ('--write-index', '-b', '-o', str(self.base / 'bad.bam'), '--order', 'store'),
                        ('--write-index', '-b', '-o', str(self.base / 'bad.bam'), '-c')):
            result = subprocess.run([self.tool, 'view', *options, str(self.view)], capture_output=True)
            self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
