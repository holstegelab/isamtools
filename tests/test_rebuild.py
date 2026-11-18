"""Exhaustive small packed-sequence fixtures for the public rebuild primitive."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PackedSequenceTests(unittest.TestCase):
    def test_all_iupac_bases_clipping_orientations_and_padding(self):
        build = Path(os.environ.get('TOOL', ROOT / 'build/isamtools')).resolve().parent
        hts = Path(os.environ.get('HTSLIB_PREFIX', build / 'vendor/htslib'))
        lmdb = Path(os.environ.get('LMDB_PREFIX', build / 'vendor/lmdb'))
        bzip = Path(os.environ.get('BZIP2_LIB', build / 'vendor/bzip2/lib/libbz2.a'))
        deflate = Path(os.environ.get('LIBDEFLATE_LIB', build / 'vendor/libdeflate/lib/libdeflate.a'))
        with tempfile.TemporaryDirectory(prefix='isam-rebuild-test-') as directory:
            executable = Path(directory) / 'rebuild_sequences'
            command = [os.environ.get('CC', 'cc'), '-O3', '-std=c11',
                       '-D_POSIX_C_SOURCE=200809L', '-no-pie', '-I' + str(ROOT / 'src'),
                       '-I' + str(hts / 'include'), str(ROOT / 'tests/rebuild_sequences.c'),
                       str(build / 'libisam.a'), str(hts / 'lib/libhts.a'),
                       str(lmdb / 'lib/liblmdb.a'), str(bzip), str(deflate),
                       '-lz', '-llzma', '-lm', '-lpthread', '-lzstd', '-o', str(executable)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            result = subprocess.run([str(executable)], check=True, capture_output=True, text=True)
            self.assertEqual(json.loads(result.stdout)['verified_cases'], 418944)


if __name__ == '__main__':
    unittest.main()
