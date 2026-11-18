#!/usr/bin/env python3
"""Run a synchronous consumer with a temporary, indexed BAM on node-local SSD.

Permanent consumer outputs must be written outside the temporary directory.
Only the unique directory created here is removed, including staged inputs.
"""
import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from reference_path import reference_path_files


def checked(command, env, capture=False):
    result = subprocess.run(command, env=env, check=True,
                            stdout=subprocess.PIPE if capture else sys.stderr)
    return result.stdout.decode().strip() if capture else None


def copy_checked(source, destination, report):
    started = time.monotonic()
    shutil.copyfile(source, destination)
    if source.stat().st_size != destination.stat().st_size:
        raise RuntimeError(f'incomplete staging copy: {source}')
    report['staging'].append({'name': destination.name, 'bytes': destination.stat().st_size,
                              'seconds': time.monotonic() - started})


def materialize(args, work, report):
    env = os.environ.copy()
    env['ISAM_SLICE_CACHE_MB'] = str(args.cache_mb)
    reference = args.reference.resolve() if args.reference else None
    reference_directory = args.reference_path.resolve() if args.reference_path else None
    bundle = args.bundle.resolve()
    cram, ssam = bundle / 'canon.cram', bundle / 'canon.ssam'
    source = cram if args.view == 'canonical' else bundle / (args.view + '.isam')
    required = [cram, source]
    if reference:
        required += [reference, Path(str(reference) + '.fai')]
    if args.view != 'canonical':
        required.append(ssam)
    if args.stage == 'all' or (args.region and args.view == 'canonical'):
        required.append(Path(str(cram) + '.crai'))
    for path in required:
        if not path.is_file():
            raise RuntimeError(f'missing required input: {path}')

    if args.view != 'canonical' and args.stage == 'none':
        for path in (ssam, source):
            if not path.resolve().is_relative_to(Path('/scratch-node')):
                raise RuntimeError(f'non-staged LMDB input {path} must resolve to /scratch-node; use --stage indexes for shared-storage indexes')

    region = [args.region] if args.region else []
    expected = args.expected_records
    inputs = work / 'inputs'
    inputs.mkdir()
    if reference_directory:
        files = reference_path_files(cram, reference_directory, args.samtools, env)
        if args.stage_reference:
            local_directory = inputs / 'reference-m5'
            local_directory.mkdir()
            for name in files:
                copy_checked(reference_directory / name, local_directory / name, report)
            reference_directory = local_directory
        env['REF_PATH'] = str(reference_directory / '%s')
        env['REF_CACHE'] = ''
        reference_options = ['--reference-path', str(reference_directory)]
    else:
        reference_options = ['-T', str(reference)]
    if args.stage_reference and reference:
        local_reference = inputs / 'reference.fa'
        copy_checked(reference, local_reference, report)
        copy_checked(Path(str(reference) + '.fai'), Path(str(local_reference) + '.fai'), report)
        reference = local_reference
        reference_options = ['-T', str(reference)]
    if args.stage != 'none' and args.view != 'canonical':
        for path in (ssam, source):
            copy_checked(path, inputs / path.name, report)
        ssam, source = inputs / ssam.name, inputs / source.name
    if args.stage == 'all':
        copy_checked(cram, inputs / cram.name, report)
        copy_checked(Path(str(cram) + '.crai'), inputs / (cram.name + '.crai'), report)
        cram = inputs / cram.name
        if args.view == 'canonical':
            source = cram

    # Open the staged LMDB files, including for counting. Cross-node readers
    # must not share LMDB's mmap-backed process-shared mutex on GPFS.
    if expected is None:
        expected = int(checked([args.isamtools, 'view', '-c'] + reference_options + [
                                '--store', str(cram), '--store-index', str(ssam), str(source)] + region,
                               env, capture=True))

    bam = work / 'input.bam'
    if args.view == 'canonical':
        producer = [args.samtools, 'view'] + (['-T', str(reference)] if reference else []) + ['-@', str(args.threads),
                    '-1', '-o', str(bam), str(source)] + region
    elif args.strategy == 'store-sort':
        producer = [args.isamtools, 'view'] + reference_options + ['-@', str(args.threads),
                    '--order', 'store', '-u', '--store', str(cram), '--store-index', str(ssam), str(source)]
    else:
        producer = [args.isamtools, 'view'] + reference_options + ['-@', str(args.threads),
                    '--order', 'coordinate', '-1', '--store', str(cram),
                    '--store-index', str(ssam), '-o', str(bam), str(source)] + region
        if args.decode_threads:
            producer[2:2] = ['--decode-threads', str(args.decode_threads)]
        if args.decode_lookahead is not None:
            producer[2:2] = ['--decode-lookahead', str(args.decode_lookahead)]
        if args.cache_selected:
            producer[2:2] = ['--cache-selected']
        if args.write_index:
            producer[2:2] = ['--write-index']
    started = time.monotonic()
    resources_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    if args.strategy == 'store-sort':
        sorter = [args.samtools, 'sort', '-@', str(args.threads), '-m', '512M', '-l', '1',
                  '-T', str(work / 'sort-tmp'), '-o', str(bam), '-']
        # Deliberately unordered output is sorted; this is not a fallback for
        # an invalid coordinate view. Check BOTH sides before indexing/consumer.
        with subprocess.Popen(producer, env=env, stdout=subprocess.PIPE) as process:
            try:
                result = subprocess.run(sorter, env=env, stdin=process.stdout, check=False)
            finally:
                process.stdout.close()
            producer_rc = process.wait()
        if producer_rc or result.returncode:
            raise RuntimeError(f'store-sort failed: producer={producer_rc}, sort={result.returncode}')
    else:
        checked(producer, env)
    report['conversion_seconds'] = time.monotonic() - started
    resources_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    report['conversion_children_user_seconds'] = resources_after.ru_utime - resources_before.ru_utime
    report['conversion_children_system_seconds'] = resources_after.ru_stime - resources_before.ru_stime
    report['conversion_children_cpu_seconds'] = (report['conversion_children_user_seconds'] +
                                                 report['conversion_children_system_seconds'])
    checked([args.samtools, 'quickcheck', '-v', str(bam)], env)
    # Indexing is also a hard coordinate-order check: do not quietly re-sort
    # an invalid migration bundle or run the consumer on an unusable BAM.
    if args.write_index:
        if not Path(str(bam) + '.bai').is_file():
            raise RuntimeError('producer did not create the requested inline BAM index')
        report['index_seconds'] = 0.0
    else:
        started = time.monotonic()
        checked([args.samtools, 'index', '-@', str(args.threads), str(bam)], env)
        report['index_seconds'] = time.monotonic() - started
    report['index_included_in_conversion'] = args.write_index
    idxstats = checked([args.samtools, 'idxstats', str(bam)], env, capture=True)
    records = sum(int(fields[2]) + int(fields[3])
                  for line in idxstats.splitlines() for fields in [line.split('\t')])
    if records != expected:
        raise RuntimeError(f'BAM record count mismatch: {records}, expected {expected}')
    report.update({'records': records, 'expected_records': expected, 'bam_bytes': bam.stat().st_size})
    return bam, env


def resolve_strategy(view, region, strategy):
    if strategy == 'auto':
        return 'coordinate' if region or view == 'canonical' else 'store-sort'
    return strategy


def execute(args):
    requested_strategy = args.strategy
    args.strategy = resolve_strategy(args.view, args.region, requested_strategy)
    if args.decode_threads and (args.strategy != 'coordinate' or args.view == 'canonical'):
        raise RuntimeError('--decode-threads requires --strategy coordinate and an iSAM view')
    if (args.cache_selected or args.decode_lookahead is not None) and not args.decode_threads:
        raise RuntimeError('--cache-selected and --decode-lookahead require --decode-threads')
    if args.write_index and (args.strategy != 'coordinate' or args.view == 'canonical'):
        raise RuntimeError('--write-index requires --strategy coordinate and an iSAM view')
    scratch = args.scratch.resolve()
    if not scratch.is_dir() or not scratch.is_relative_to(Path('/scratch-node')):
        raise RuntimeError('--scratch must be an existing node-local /scratch-node directory')
    work = Path(tempfile.mkdtemp(prefix='isam-with-bam-', dir=scratch))
    report = {'bundle': str(args.bundle.resolve()), 'view': args.view,
              'reference_mode': 'm5' if args.reference_path else 'fasta',
              'stage': args.stage, 'requested_strategy': requested_strategy,
              'strategy': args.strategy, 'decode_threads': args.decode_threads,
              'decode_lookahead': args.decode_lookahead, 'cache_selected': args.cache_selected,
              'staging': [], 'scratch': str(work)}
    try:
        bam, env = materialize(args, work, report)
        command = [part.replace('{bam}', str(bam)).replace('{scratch}', str(work)) for part in args.command]
        started = time.monotonic()
        consumer = subprocess.run(command, env=env, check=False)
        report['consumer_seconds'] = time.monotonic() - started
        report['consumer_returncode'] = consumer.returncode
        return consumer.returncode if consumer.returncode >= 0 else 128 - consumer.returncode
    except BaseException as error:
        report['error'] = str(error)
        raise
    finally:
        shutil.rmtree(work)
        report['local_bulk_data_removed'] = True
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps(report, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', required=True, type=Path)
    parser.add_argument('--view', required=True, help='View name without .isam, or canonical')
    references = parser.add_mutually_exclusive_group(required=True)
    references.add_argument('--reference', type=Path, help='Canonical hg38 FASTA, also for CHM13 views')
    references.add_argument('--reference-path', type=Path,
                            help='Verified raw M5 directory from reference_path.py; requires a current reader')
    parser.add_argument('--scratch', required=True, type=Path)
    parser.add_argument('--isamtools', default='isamtools')
    parser.add_argument('--samtools', default='samtools')
    parser.add_argument('--threads', type=int, default=3)
    parser.add_argument('--decode-threads', type=int, default=0,
                        help='Opt-in independent CRAM container workers (1..64); requires --strategy coordinate')
    parser.add_argument('--decode-lookahead', type=int,
                        help='Coordinate metadata window (1..65536); requires --decode-threads')
    parser.add_argument('--cache-selected', action='store_true',
                        help='Plan requested payloads and release final uses; requires --decode-threads')
    parser.add_argument('--write-index', action='store_true',
                        help='Build BAI during coordinate iSAM export, avoiding a separate BAM scan')
    parser.add_argument('--cache-mb', type=int, default=4096)
    parser.add_argument('--stage', choices=['none', 'indexes', 'all'], default='indexes',
                        help='Default: private SSD copies of LMDB indexes; none is for local inputs only')
    parser.add_argument('--stage-reference', action='store_true')
    parser.add_argument('--strategy', choices=['auto', 'coordinate', 'store-sort'], default='auto',
                        help='auto: store-sort for whole iSAM views, coordinate for regions/canonical; coordinate: indexed direct output; store-sort: whole canonical stream and SSD sort')
    parser.add_argument('--region')
    parser.add_argument('--expected-records', type=int, help='Alignment record count for the selected scope (whole file or --region), not unique reads; required for canonical')
    parser.add_argument('--report', type=Path)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ['--']:
        args.command = args.command[1:]
    if not args.command or not any('{bam}' in part for part in args.command):
        parser.error('provide a synchronous consumer command after -- containing {bam}')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', args.view):
        parser.error('invalid view name')
    if args.threads < 1 or args.cache_mb < 1:
        parser.error('threads and cache-mb must be positive')
    if args.decode_threads < 0 or args.decode_threads > 64:
        parser.error('decode-threads must be 0 (disabled) or 1..64')
    if args.decode_lookahead is not None and not 1 <= args.decode_lookahead <= 65536:
        parser.error('decode-lookahead must be 1..65536')
    if args.expected_records is not None and args.expected_records < 0:
        parser.error('expected-records cannot be negative')
    if args.view == 'canonical' and args.expected_records is None:
        parser.error('canonical requires --expected-records from its validated index/migration receipt')
    if args.strategy == 'store-sort' and (args.region or args.view == 'canonical'):
        parser.error('store-sort requires a whole iSAM view (no --region and not canonical)')
    if args.view == 'canonical' and args.stage == 'indexes':
        # Canonical CRAM has no LMDB input, so the safe default needs no staging.
        args.stage = 'none'
    if args.report and args.report.exists():
        parser.error('report path already exists; refusing overwrite')
    try:
        return execute(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f'materialization failed: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
