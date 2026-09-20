"""Approximate exon-query mapping used only to seed local BLAST windows."""
import os
import subprocess
import sys


def iter_minimap_candidates(args, work_dir):
    index = os.path.join(work_dir, 'assembly.k19.w5.mmi')
    command = [args.minimap2, '-k19', '-w5', '-I8G', '--idx-no-seq',
               '-t', str(args.threads), '-d', index, args.query]
    print('Indexing minimap2 candidates:', ' '.join(command), file=sys.stderr)
    subprocess.run(command, check=True)
    command = [args.minimap2, '-P', '--secondary=yes', '-n3', '-m40',
               '-f1000', '--q-occ-frac=0', '--no-long-join', '-g1000',
               '-K', str(args.minimap_batch_bases), '-t', str(args.threads),
               index, args.exon_fasta or f'{args.db}.exons.fa']
    print('Minimap2 candidate discovery (approximate PAF; no identity/MAPQ filter):',
          ' '.join(command), file=sys.stderr)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    try:
        yield from proc.stdout
        code = proc.wait()
        if code:
            raise subprocess.CalledProcessError(code, command)
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.terminate()
            proc.wait()
        if os.path.exists(index):
            os.unlink(index)


def parse_candidate(line):
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 12:
        raise ValueError('Malformed minimap2 PAF record')
    qid, contig = fields[0], fields[5]
    length, start, end = int(fields[6]), int(fields[7]), int(fields[8])
    if not 0 <= start < end <= length:
        raise ValueError(f'Invalid minimap2 target coordinates: {contig}:{start}-{end}/{length}')
    return qid, contig, length, start, end
