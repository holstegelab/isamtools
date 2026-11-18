#!/usr/bin/env python3
"""Prepare and inspect verified raw M5 references for local HTSlib REF_PATH.

Preparation requires pysam; inspection and the BAM consumer helper use only
the standard library. Input CRAM and FASTA files are never modified.
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def require_equal(actual, expected, label):
    if actual != expected:
        raise RuntimeError(f'{label}: actual={actual}, expected={expected}')


def read_reference_sequences(cram, samtools, env=None):
    result = subprocess.run([samtools, 'view', '-H', str(cram)], check=True,
                            capture_output=True, text=True, env=env)
    sequences = []
    for line in result.stdout.splitlines():
        if not line.startswith('@SQ\t'):
            continue
        fields = dict(field.split(':', 1) for field in line.split('\t')[1:])
        if not fields.get('SN') or not fields.get('LN', '').isdigit() or \
                not re.fullmatch('[0-9a-fA-F]{32}', fields.get('M5', '')):
            raise RuntimeError('canonical SQ requires SN, positive LN and M5')
        length = int(fields['LN'])
        if length <= 0:
            raise RuntimeError('canonical SQ length must be positive')
        sequences.append({'SN': fields['SN'], 'LN': length, 'M5': fields['M5']})
    if not sequences:
        raise RuntimeError('canonical CRAM header has no SQ references')
    return sequences


def reference_path_files(cram, directory, samtools, env=None):
    directory = directory.resolve()
    # These characters have special REF_PATH list/substitution semantics.
    if not directory.is_dir() or ':' in str(directory) or '%' in str(directory):
        raise RuntimeError('reference-path must be an existing directory without : or %')
    files = {}
    for sequence in read_reference_sequences(cram, samtools, env):
        checksum, length = sequence['M5'], sequence['LN']
        path = directory / checksum
        if not path.is_file() or path.stat().st_size != length:
            raise RuntimeError('missing/wrong-length M5 reference: ' + str(path))
        if checksum in files:
            require_equal(files[checksum], length, 'duplicate M5 length')
        files[checksum] = length
    return files


def populate_reference_path(reference, sequences, destination, chunk_bases=8 * 1024 * 1024):
    """Create fresh raw M5 files, verifying every canonical contig in chunks."""
    import pysam

    if chunk_bases < 1 or not sequences:
        raise RuntimeError('positive reference chunk size and nonempty SQ header required')
    if ':' in str(destination) or '%' in str(destination):
        raise RuntimeError('reference-path directory cannot contain : or %')
    # Exclusive creation is outside cleanup: an existing directory belongs to
    # the caller, and must not be removed if this preparation refuses it.
    destination.mkdir()
    checksums = {}
    try:
        with pysam.FastaFile(str(reference)) as fasta:
            for sequence in sequences:
                name, length = sequence['SN'], sequence['LN']
                # HTSlib substitutes the literal header M5. Preserve filename
                # case; only the digest comparison is case-insensitive.
                checksum = sequence.get('M5', '')
                if length <= 0 or not re.fullmatch('[0-9a-fA-F]{32}', checksum):
                    raise RuntimeError('invalid or missing canonical SQ length/M5: ' + name)
                require_equal(fasta.get_reference_length(name), length, name + ' reference length')
                temporary = destination / (checksum + '.tmp')
                digest = hashlib.md5()
                with temporary.open('xb') as output:
                    for start in range(0, length, chunk_bases):
                        end = min(length, start + chunk_bases)
                        bases = fasta.fetch(name, start, end).encode('ascii').upper()
                        require_equal(len(bases), end - start, name + ' reference chunk length')
                        digest.update(bases)
                        output.write(bases)
                require_equal(digest.hexdigest(), checksum.lower(), name + ' canonical reference M5')
                if checksum in checksums:
                    require_equal(length, checksums[checksum], name + ' duplicate M5 length')
                    temporary.unlink()
                else:
                    temporary.replace(destination / checksum)
                    checksums[checksum] = length
    except BaseException:
        shutil.rmtree(destination)
        raise
    return {'contigs_verified': len(sequences), 'unique_m5_files': len(checksums),
            'bytes': sum(checksums.values()), 'chunk_bases': chunk_bases,
            'verification': 'every canonical SQ length and M5 versus uppercase FASTA bases'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cram', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path, help='Fresh directory; existing paths are refused')
    parser.add_argument('--samtools', default='samtools')
    args = parser.parse_args()
    try:
        sequences = read_reference_sequences(args.cram, args.samtools)
        report = populate_reference_path(args.reference, sequences, args.output)
        print(json.dumps(dict(report, directory=str(args.output.resolve()))))
        return 0
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f'reference preparation failed: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
