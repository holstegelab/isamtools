"""Native raw-reference preparation, validation and one-step BAM conversion."""
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ReferenceIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tool = str(Path(os.environ.get('TOOL', ROOT / 'build/isamtools')).resolve())
        cls.samtools = shutil.which('samtools')
        if not cls.samtools or not Path(cls.tool).is_file():
            raise RuntimeError('samtools and a built isamtools are required')
        cls.fixture = tempfile.TemporaryDirectory(prefix='isam-native-reference-')
        cls.base = Path(cls.fixture.name)
        cls.ref = cls.base / 'reference.fa'
        cls.sequences = {'chr1': 'acgt' * 64, 'copy': 'acgt' * 64, 'chr2': 'tgca' * 64}
        cls.ref.write_text(''.join('>' + n + '\n' + s + '\n' for n, s in cls.sequences.items()))
        header = '@HD\tSO:coordinate\n' + ''.join(f'@SQ\tSN:{n}\tLN:{len(s)}\n' for n, s in cls.sequences.items())
        header += '@RG\tID:movie\tSM:fixture\tPL:PACBIO\n'
        records = [header]
        for i in range(8):
            records.append(f'read{i}\t0\tchr1\t{1+4*i}\t60\t32M\t*\t0\t0\t' +
                           'ACGT' * 8 + '\t' + 'I' * 32 + '\tRG:Z:movie\tHP:i:2\tPS:i:1\tMM:Z:C+m,0;\tML:B:C,200\n')
        records.append('unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t' + 'ACGT' * 8 + '\t' + 'I' * 32 + '\tRG:Z:movie\tHP:i:1\n')
        sam = cls.base / 'source.sam'
        sam.write_text(''.join(records))
        cls.bam, cls.cram, cls.ssam, cls.view = [cls.base / n for n in ('source.bam', 'canon.cram', 'canon.ssam', 'view.isam')]
        cls.run_ok([cls.samtools, 'faidx', str(cls.ref)])
        cls.run_ok([cls.samtools, 'view', '-b', '-o', str(cls.bam), str(sam)])
        cls.run_ok([cls.samtools, 'view', '-C', '-T', str(cls.ref), '-o', str(cls.cram), str(cls.bam)])
        cls.run_ok([cls.samtools, 'index', str(cls.cram)])
        cls.run_ok([cls.tool, 'index', '-T', str(cls.ref), '-o', str(cls.ssam), str(cls.cram)])
        cls.run_ok([cls.tool, 'import', '--fast-verify', '--store', str(cls.cram), '--store-index', str(cls.ssam),
                    '--aln', str(cls.view), str(cls.bam)])
        cls.run_ok([cls.tool, 'sort', '--all', str(cls.view)])

    @classmethod
    def tearDownClass(cls):
        cls.fixture.cleanup()

    @staticmethod
    def run_ok(cmd):
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if r.returncode:
            raise RuntimeError(f'{cmd}: {r.stderr}')
        return r

    def prepare(self, output, ref=None, cram=None):
        return [self.tool, 'ref-index', '--cram', str(cram or self.cram), '-T', str(ref or self.ref), '-o', str(output)]

    def test_native_command_checks_and_deduplicates_every_contig(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            out = Path(d) / 'raw'
            r = self.run_ok(self.prepare(out))
            expected = {hashlib.md5(s.upper().encode()).hexdigest(): s.upper().encode() for s in self.sequences.values()}
            self.assertEqual({p.name for p in out.iterdir()}, set(expected))
            for name, sequence in expected.items():
                self.assertEqual((out / name).read_bytes(), sequence)
            self.assertIn('verified 3 contigs, wrote 2 files', r.stderr)

    def test_existing_directory_and_symlink_are_preserved(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            out = Path(d) / 'owned'
            out.mkdir()
            marker = out / 'mine'
            marker.write_text('keep')
            link = Path(d) / 'link'
            link.symlink_to(out, target_is_directory=True)
            for path in (out, link):
                r = subprocess.run(self.prepare(path), capture_output=True, text=True)
                self.assertNotEqual(r.returncode, 0)
                self.assertEqual(marker.read_text(), 'keep')
                self.assertEqual(list(out.iterdir()), [marker])

    def test_checksum_failure_removes_only_own_outputs(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            wrong = Path(d) / 'wrong.fa'
            wrong.write_text(''.join('>' + n + '\n' + (('A' + s[1:]) if n == 'chr2' else s) + '\n'
                                     for n, s in self.sequences.items()))
            out = Path(d) / 'failed'
            r = subprocess.run(self.prepare(out, ref=wrong), capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn('MD5 differs', r.stderr)
            self.assertFalse(out.exists())
            self.assertTrue(wrong.exists())

    def test_native_preparation_creates_missing_fai(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            ref = Path(d) / 'new.fa'
            shutil.copyfile(self.ref, ref)
            self.assertFalse(Path(str(ref) + '.fai').exists())
            self.run_ok(self.prepare(Path(d) / 'raw', ref=ref))
            self.assertTrue(Path(str(ref) + '.fai').exists())

    def test_chunk_boundary_has_no_missing_or_duplicated_bases(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            root = Path(d)
            sequence = 'acgt' * (2 * 1024 * 1024) + 'acgtn'
            ref = root / 'large.fa'
            ref.write_text('>large\n' + sequence + '\n')
            sam = root / 'large.sam'
            sam.write_text(f'@SQ\tSN:large\tLN:{len(sequence)}\nread\t0\tlarge\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n')
            cram = root / 'large.cram'
            self.run_ok([self.samtools, 'view', '-C', '-T', str(ref), '-o', str(cram), str(sam)])
            out = root / 'raw'
            self.run_ok(self.prepare(out, ref=ref, cram=cram))
            checksum = hashlib.md5(sequence.upper().encode()).hexdigest()
            self.assertEqual((out / checksum).read_bytes(), sequence.upper().encode())

    def test_invalid_m5_is_reported_and_preparation_is_cleaned(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            root = Path(d)
            header = self.run_ok([self.samtools, 'view', '-H', str(self.cram)]).stdout
            header = re.sub(r'\tM5:[0-9a-fA-F]{32}', '\tM5:invalid', header)
            h = root / 'header.sam'
            h.write_text(header)
            cram = root / 'missing.cram'
            with cram.open('wb') as output:
                subprocess.run([self.samtools, 'reheader', str(h), str(self.cram)], stdout=output, check=True)
            out = root / 'raw'
            result = subprocess.run(self.prepare(out, cram=cram), capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('M5 checksum', result.stderr)
            self.assertFalse(out.exists())

    def test_literal_uppercase_m5_filenames(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            header = self.run_ok([self.samtools, 'view', '-H', str(self.cram)]).stdout
            header = re.sub(r'M5:([0-9a-f]{32})', lambda m:'M5:' + m.group(1).upper(), header)
            h = Path(d) / 'header.sam'
            h.write_text(header)
            cram = Path(d) / 'upper.cram'
            with cram.open('wb') as output:
                subprocess.run([self.samtools, 'reheader', str(h), str(self.cram)], stdout=output, check=True)
            out = Path(d) / 'raw'
            self.run_ok(self.prepare(out, cram=cram))
            self.assertEqual({p.name for p in out.iterdir()},
                             {hashlib.md5(s.upper().encode()).hexdigest().upper() for s in self.sequences.values()})

    def test_one_step_conversion_preserves_original_records_and_tags(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            out = Path(d) / 'raw'
            prefix = [self.tool, 'view', '-h', '--store', str(self.cram), '--store-index', str(self.ssam)]
            baseline = self.run_ok(prefix + ['-T', str(self.ref), str(self.view)])
            self.assertIn('ref-index', baseline.stderr)
            prepared = self.run_ok(prefix + ['-T', str(self.ref), '--prepare-reference', str(out),
                                             '--decode-threads', '2', str(self.view)])
            self.assertEqual(prepared.stdout, baseline.stdout)
            self.assertIn('verified 3 contigs', prepared.stderr)
            reused = self.run_ok(prefix + ['--reference-path', str(out), str(self.view)])
            self.assertEqual(reused.stdout, baseline.stdout)
            self.assertNotIn('FASTA coordinate decoding', reused.stderr)
            self.assertNotIn('FASTA coordinate decoding', prepared.stderr)

    def test_preparation_options_and_metadata_queries_have_no_side_effects(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            out = Path(d) / 'raw'
            for opts in ([], ['-T', str(self.ref), '-c'], ['-T', str(self.ref), '-H'],
                         ['-T', str(self.ref), '-C'], ['-T', str(self.ref), '--reference-path', str(out)]):
                r = subprocess.run([self.tool, 'view', '--prepare-reference', str(out)] + opts + [str(self.view)],
                                   capture_output=True, text=True)
                self.assertNotEqual(r.returncode, 0)
                self.assertFalse(out.exists())
            count = self.run_ok([self.tool, 'view', '-c', str(self.view)])
            self.assertEqual(count.stdout.strip(), '9')
            self.assertNotIn('FASTA coordinate decoding', count.stderr)

    def test_preparation_rejects_wrong_length_and_special_directory(self):
        with tempfile.TemporaryDirectory(dir=self.base) as d:
            short = Path(d) / 'short.fa'
            short.write_text('>chr1\nACGT\n')
            out = Path(d) / 'raw'
            r = subprocess.run(self.prepare(out, ref=short), capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn('length differs', r.stderr)
            self.assertFalse(out.exists())
            for name in ('bad:directory', 'bad%directory'):
                path = Path(d) / name
                r = subprocess.run(self.prepare(path), capture_output=True, text=True)
                self.assertNotEqual(r.returncode, 0)
                self.assertFalse(path.exists())


if __name__ == '__main__':
    unittest.main()
