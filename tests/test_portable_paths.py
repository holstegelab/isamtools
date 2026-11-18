import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class PortableBundleTests(unittest.TestCase):
    def setUp(self):
        self.samtools = shutil.which('samtools')
        project = Path(__file__).resolve().parents[1]
        self.tool = Path(os.environ.get('TOOL', project / 'build/isamtools')).resolve()
        if self.samtools is None or not self.tool.is_file():
            self.fail('samtools and a built isamtools are required')
        self.fixture = tempfile.TemporaryDirectory(prefix='isam-portable-bundle-test-')
        self.addCleanup(self.fixture.cleanup)
        self.base = base = Path(self.fixture.name)
        self.bundle = bundle = base / 'old_worker_scratch'
        bundle.mkdir()
        self.other = base / 'unrelated_cwd'
        self.other.mkdir()
        self.reference = reference = base / 'reference.fa'
        reference.write_text('>chr1\n' + 'ACGT' * 250 + '\n')
        source = base / 'source.sam'
        source.write_text(
                '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:movie\tSM:synthetic\n'
                'read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:movie\tHP:i:1\tPS:i:17\tMM:Z:C+m,0;\tML:B:C,200\tzm:i:5\n'
                'read1\t2048\tchr1\t11\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:movie\tHP:i:2\tPS:i:23\tMM:Z:C+m,0;\tML:B:C,200\tzm:i:5\n'
                'read2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:movie\tMM:Z:C+m,0;\tML:B:C,201\tzm:i:6\n')
        self.bam = bam = base / 'source.bam'
        self.view = view = bundle / 'hifi.hg38.isam'
        tool = str(self.tool)
        for command in (
            [self.samtools, 'faidx', str(reference)],
            [self.samtools, 'view', '-b', '-o', str(bam), str(source)],
            [self.samtools, 'view', '-C', '-T', str(reference), '-o', 'old_worker_scratch/canon.cram', str(source)],
            [self.samtools, 'index', 'old_worker_scratch/canon.cram'],
            [tool, 'index', '-T', str(reference), 'old_worker_scratch/canon.cram'],
            [tool, 'import', '--strict-tag-policy', '--fast-verify', '--store', 'old_worker_scratch/canon.cram',
             '--store-index', 'old_worker_scratch/canon.ssam', '--aln', str(view), str(bam)],
            [tool, 'sort', '--all', str(view)],
        ):
            self.run_command(command, cwd=base)

    def run_command(self, command, cwd=None):
        result = subprocess.run(command, text=True, capture_output=True, cwd=cwd)
        self.assertEqual(result.returncode, 0, f'{command!r}\n{result.stderr}')
        return result.stdout

    def verify_reconstruction(self, view, extra=(), cwd=None):
        restored = self.base / 'restored.bam'
        self.run_command([str(self.tool), 'view', '-T', str(self.reference), '-1', '-o', str(restored),
                          *extra, str(view)], cwd=cwd)
        self.run_command([self.samtools, 'quickcheck', '-v', str(restored)])
        self.run_command([self.samtools, 'index', str(restored)])
        counts = self.run_command([self.samtools, 'idxstats', str(restored)])
        self.assertEqual(sum(int(fields[2]) + int(fields[3])
                             for line in counts.splitlines() for fields in [line.split('\t')]), 3)

        def records(text):
            return [(tuple(fields[:11]), tuple(sorted(fields[11:])))
                    for line in text.splitlines() for fields in [line.split('\t')]]

        original = self.run_command([self.samtools, 'view', str(self.bam)])
        observed = self.run_command([self.samtools, 'view', str(restored)])
        self.assertEqual(records(observed), records(original))

    def move_bundle(self):
        self.run_command([str(self.tool), 'update', '--store-path', 'canon.cram',
                          '--store-index-path', 'canon.ssam', str(self.view)])
        relocated = self.base / 'restored_archive_bundle'
        self.bundle.rename(relocated)
        self.assertFalse(self.bundle.exists())
        return relocated / 'hifi.hg38.isam'

    def test_cwd_relative_import_paths_are_persisted_without_rebasing(self):
        self.verify_reconstruction(self.view, cwd=self.other)

    def test_relative_archive_paths_after_old_scratch_is_gone(self):
        view = self.move_bundle()
        self.verify_reconstruction(view, cwd=self.other)
        self.verify_reconstruction(Path('restored_archive_bundle/hifi.hg38.isam'), cwd=self.base)

    def test_explicit_relative_overrides_still_use_caller_cwd(self):
        view = self.move_bundle()
        self.run_command([str(self.tool), 'update', '--store-path', 'absent.cram',
                          '--store-index-path', 'absent.ssam', str(view)])
        self.verify_reconstruction(view, extra=('--store', 'restored_archive_bundle/canon.cram',
                                               '--store-index', 'restored_archive_bundle/canon.ssam'), cwd=self.base)

    def test_query_resolves_relative_store_path_from_relocated_ssam(self):
        self.run_command([str(self.tool), 'update', '--set', 'store_path=canon.cram',
                          str(self.bundle / 'canon.ssam')])
        view = self.move_bundle()
        names = self.base / 'qnames.txt'
        names.write_text('read1\nread2\n')
        observed = self.run_command([str(self.tool), 'query', '--source',
                                    'canon:' + str(view.parent / 'canon.ssam'),
                                    '--qnames', str(names), '--columns',
                                    'qname,canon.region', '--no-header'], cwd=self.other)
        self.assertEqual(observed.splitlines(), ['read1\tchr1:1-8', 'read2\t.'])

    def test_rehint_uses_bundle_paths_and_persists_cwd_relative_overrides(self):
        view = self.move_bundle()
        self.run_command([str(self.tool), 'rehint', str(view)], cwd=self.other)
        self.verify_reconstruction(view, cwd=self.other)
        twice = self.base / 'twice_restored_bundle'
        view.parent.rename(twice)
        view = twice / 'hifi.hg38.isam'
        self.verify_reconstruction(view, cwd=self.other)
        self.run_command([str(self.tool), 'rehint', '--store', 'twice_restored_bundle/canon.cram',
                          '--store-index', 'twice_restored_bundle/canon.ssam', str(view)], cwd=self.base)
        self.verify_reconstruction(view, cwd=self.other)


if __name__ == '__main__':
    unittest.main()
