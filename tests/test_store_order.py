"""Sequential canonical materialization must equal the indexed view path."""
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class StoreOrderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tool = str(Path(os.environ.get('TOOL', ROOT / 'build/isamtools')).resolve())
        cls.samtools = shutil.which('samtools')
        if not cls.samtools or not Path(cls.tool).is_file():
            raise RuntimeError('samtools and built isamtools are required')
        cls.fixture = tempfile.TemporaryDirectory(prefix='isam-store-order-')
        cls.base = Path(cls.fixture.name)
        cls.reference = cls.base / 'reference.fa'
        cls.reference.write_text('>chr1\n' + 'ACGT' * 250 + '\n')
        header = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n' \
                 '@RG\tID:canon\tSM:synthetic\n@RG\tID:view\tSM:synthetic\n'

        def record(name, flag, pos, rg='canon', hp=1, cigar='8M', seq='ACGTTGCA', qual='IHGFEDCB'):
            return '\t'.join([name, str(flag), 'chr1' if pos else '*', str(pos), '60',
                              cigar if pos else '*', '*', '0', '0', seq, qual,
                              'RG:Z:' + rg, 'HP:i:' + str(hp), 'PS:i:17',
                              'MM:Z:C+m,0;', 'ML:B:C,200', 'zm:i:101']) + '\n'

        # qA upgrades an earlier supplementary locator to its later primary;
        # qSupp is supplementary-only, with two potential full-length payloads.
        canonical = header + record('qA', 2048, 10, hp=2) + record('qSupp', 2048, 20) \
                    + record('qSupp', 2048, 30, hp=2) \
                    + record('qHard', 2048, 40, cigar='2H6M', seq='GTTGCA', qual='GFEDCB') \
                    + record('qA', 0, 80) + record('qA', 256, 90) + record('qUnmap', 4, 0)
        alignment = header + record('qA', 16, 5, 'view', 2, seq='TGCAACGT', qual='BCDEFGHI') \
                    + record('qA', 2048, 15, 'view', 2) + record('qSupp', 0, 25, 'view') \
                    + record('qSupp', 256, 35, 'view', 2) \
                    + record('qA', 2048, 45, 'view', 2, '2H6M', 'GTTGCA', 'GFEDCB') \
                    + record('qUnmap', 4, 0, 'view', 2)
        (cls.base / 'canonical.sam').write_text(canonical)
        (cls.base / 'alignment.sam').write_text(alignment)
        (cls.base / 'store_tags.tsv').write_text('qA\tZS:Z:store_annotation\tZX:i:31\n')
        (cls.base / 'view_tags.tsv').write_text('qA\tZV:Z:view_annotation\tZX:i:45\n')
        cls.checked([cls.samtools, 'faidx', str(cls.reference)])
        cls.checked([cls.samtools, 'view', '-b', '-o', str(cls.base / 'alignment.bam'),
                 str(cls.base / 'alignment.sam')])
        cls.views = {}
        for fmt in ('bam', 'cram'):
            store = cls.base / ('canonical.' + fmt)
            view = cls.base / (fmt + '.isam')
            cls.checked([cls.samtools, 'view', '-b' if fmt == 'bam' else '-C',
                     '-T', str(cls.reference), '-o', str(store), str(cls.base / 'canonical.sam')])
            if fmt == 'cram':
                cls.checked([cls.samtools, 'index', str(store)])
            index = cls.base / (fmt + '.ssam')
            cls.checked([cls.tool, 'index', '-T', str(cls.reference), '-o', str(index), str(store)])
            cls.checked([cls.tool, 'import', '--fast-verify', '--store', str(store),
                     '--store-index', str(index), '--aln', str(view), str(cls.base / 'alignment.bam')])
            cls.checked([cls.tool, 'sort', '--all', str(view)])
            cls.checked([cls.tool, 'annotate', str(index), '--from', str(cls.base / 'store_tags.tsv')])
            cls.checked([cls.tool, 'annotate', str(view), '--from', str(cls.base / 'view_tags.tsv')])
            cls.views[fmt] = (store, index, view)

    @classmethod
    def tearDownClass(cls):
        cls.fixture.cleanup()

    @staticmethod
    def checked(command):
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode:
            raise RuntimeError(f'{command} failed: {result.stderr}')
        return result.stdout

    def view(self, fmt, order, threads=3, options=()):
        store, index, view = self.views[fmt]
        return self.checked([self.tool, 'view', '-h', '-T', str(self.reference), '-@', str(threads),
                         '--order', order, '--store', str(store), '--store-index', str(index),
                         *options, str(view)])

    @staticmethod
    def records(text):
        return sorted(line for line in text.splitlines() if not line.startswith('@'))

    def test_bam_and_cram_exact_payloads_tags_and_header(self):
        for fmt in self.views:
            expected = self.records(self.view(fmt, 'coordinate'))
            self.assertEqual(len(expected), 6)
            for threads in (1, 4):
                with self.subTest(format=fmt, threads=threads):
                    output = self.view(fmt, 'store', threads)
                    self.assertEqual(self.records(output), expected)
                    self.assertIn('SO:unsorted', output.splitlines()[0])
                    for record in self.records(output):
                        for tag in ('MM:Z:C+m,0;', 'ML:B:C,200', 'PS:i:17', 'RG:Z:view', 'zm:i:101'):
                            self.assertIn(tag, record)
                        if record.startswith('qA\t'):
                            for tag in ('ZS:Z:store_annotation', 'ZV:Z:view_annotation', 'ZX:i:45'):
                                self.assertIn(tag, record)
                            self.assertNotIn('ZX:i:31', record)

    def test_filters_are_unchanged(self):
        for fmt in self.views:
            for options in (('-f', '4'), ('-F', '256'), ('-r', 'view'), ('-d', 'HP:2'),
                            ('-p', '-f', '4'), ('-x', 'MM,ML')):
                with self.subTest(format=fmt, options=options):
                    self.assertEqual(self.records(self.view(fmt, 'store', options=options)),
                                     self.records(self.view(fmt, 'coordinate', options=options)))

    def test_region_and_unknown_order_rejected(self):
        _, _, view = self.views['cram']
        for options in (('--order', 'store', str(view), 'chr1:1-100'),
                        ('--order', 'invalid', str(view))):
            result = subprocess.run([self.tool, 'view', *options], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)

    def test_stream_then_sort_is_indexable(self):
        store, index, view = self.views['cram']
        raw, ordered = self.base / 'stream.bam', self.base / 'sorted.bam'
        self.checked([self.tool, 'view', '--order', 'store', '-@', '4', '-u', '-T', str(self.reference),
                  '--store', str(store), '--store-index', str(index), '-o', str(raw), str(view)])
        self.checked([self.samtools, 'sort', '-o', str(ordered), str(raw)])
        self.checked([self.samtools, 'index', str(ordered)])
        self.assertEqual(self.checked([self.samtools, 'view', '-c', str(ordered)]).strip(), '6')

    def test_coordinate_decoders_keep_exact_locators_clipping_and_annotations(self):
        expected = self.view('cram', 'coordinate')
        for workers in (1, 2, 4):
            with self.subTest(workers=workers):
                self.assertEqual(self.view('cram', 'coordinate',
                                           options=('--decode-threads', str(workers))), expected)


if __name__ == '__main__':
    unittest.main()
