"""Disk-backed coarse-to-local exon searches, with per-gene target windows."""
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import closing
from dataclasses import dataclass, field
import os
import pickle
import sqlite3
import subprocess
import sys
import tempfile
import threading

from minimap_candidates import parse_candidate


class LocalProcesses:
    """Track only this local-search pool's children, including extraction jobs."""
    def __init__(self):
        self.lock = threading.Lock()
        self.children = set()
        self.stopped = False

    def start(self, command, **kwargs):
        with self.lock:
            if self.stopped:
                raise RuntimeError('Local search cancelled after another gene failed')
            proc = subprocess.Popen(command, **kwargs)
            self.children.add(proc)
            return proc

    def release(self, proc):
        with self.lock:
            self.children.discard(proc)

    def run(self, command, **kwargs):
        proc = self.start(command, **kwargs)
        try:
            code = proc.wait()
            if code:
                raise subprocess.CalledProcessError(code, command)
        finally:
            self.release(proc)

    def cancel(self):
        with self.lock:
            self.stopped = True
            for proc in self.children:
                if proc.poll() is None:
                    try:
                        proc.terminate()
                    except ProcessLookupError:
                        pass


def gene_key(alias):
    names = tuple(sorted(x for x in (v.strip() for v in alias.gene_name.split(',')) if x and x != '.'))
    ids = tuple(sorted(x for x in (v.strip() for v in alias.gene_id.split(',')) if x and x != '.'))
    return ('name', names, alias.chrom, alias.strand) if names else ('id', ids, alias.chrom, alias.strand)


@dataclass
class Gene:
    start: int
    end: int
    queries: dict = field(default_factory=dict)
    exons: set = field(default_factory=set)

    @property
    def padding(self):
        # ceil(1.5 * genomic span), on EACH side of every hit.
        return (3 * (self.end - self.start) + 1) // 2


def collect_genes(alias_map):
    genes, keys, query_genes = [], {}, defaultdict(set)
    seen_lists = set()
    for qid, aliases in alias_map.items():
        if id(aliases) in seen_lists:
            continue
        seen_lists.add(id(aliases))
        for alias in aliases:
            key = gene_key(alias)
            if not key[1]:
                raise ValueError(f'Exon {alias.exon_id_full} has neither gene name nor gene ID')
            start = alias.gene_start0 if alias.gene_start0 is not None else alias.start0
            end = alias.gene_end0 if alias.gene_end0 is not None else alias.end0
            if key not in keys:
                keys[key] = len(genes)
                genes.append(Gene(start, end))
            gid = keys[key]
            gene = genes[gid]
            gene.start, gene.end = min(gene.start, start), max(gene.end, end)
            gene.queries.setdefault(qid, []).append(alias)
            gene.exons.add(alias.exon_id_full)
            query_genes[qid].add(gid)
    return genes, query_genes


def merged_intervals(intervals):
    """Merge a sorted stream of zero-based intervals; touching windows join."""
    start = end = None
    for left, right in intervals:
        if start is None:
            start, end = left, right
        elif left <= end:
            end = max(end, right)
        else:
            yield start, end
            start, end = left, right
    if start is not None:
        yield start, end


def balanced_hits(expected, hits):
    """Count distinct overlapping hit clusters per original exon and strand.

    hits must be sorted by exon, strand, start, end. Duplicate or overlapping
    HSPs at a single physical copy count once; adjacent copies remain separate.
    """
    counts = dict.fromkeys(expected, 0)
    previous = None
    end = -1
    for exon, strand, left, right in hits:
        key = (exon, strand)
        if key != previous or left >= end:
            counts[exon] += 1
            end = right
        else:
            end = max(end, right)
        previous = key
    values = set(counts.values())
    return len(values) == 1 and next(iter(values), 0) > 0


def passes_filters(row, args):
    coverage = args.min_exon_coverage
    if args.merged_exon_queries and args.qcov_hsp_perc is not None:
        coverage = max(coverage, args.qcov_hsp_perc)
    return (row.selection_AS > args.min_as
            and row.selection_percent_identity > args.min_identity
            and (not args.merged_exon_queries or row.selection_percent_identity >= args.blast_perc_identity)
            and row.selection_coverage >= coverage)


def index_queries(api, source, path):
    """Normalize once on disk, retaining masking, for random access by gene."""
    offsets = {}
    with open(path, 'wb') as out:
        for header, sequence in api.iter_query_records(source):
            qid = header[1:].split()[0]
            record = f'>{qid}\n{sequence}\n'.encode()
            offsets[qid] = (out.tell(), len(record))
            out.write(record)
    return offsets


def extract_windows(args, api, assembly_db, work_dir, windows, processes=None):
    """Extract forward target windows in request order and give them stable IDs."""
    requests = os.path.join(work_dir, 'window_requests.txt')
    extracted = os.path.join(work_dir, 'extracted_windows.fa')
    subject = os.path.join(work_dir, 'local_windows.fa')
    with open(requests, 'w') as out:
        for contig, start, end in windows:
            out.write(f'{contig} {start + 1}-{end}\n')
    command = [args.blastdbcmd, '-db', assembly_db, '-entry_batch', requests,
               '-outfmt', '%f', '-out', extracted]
    if processes is None:
        subprocess.run(command, check=True)
    else:
        processes.run(command)
    count = 0
    with open(subject, 'w') as out:
        for i, (_, sequence) in enumerate(api.iter_query_records(extracted)):
            if i >= len(windows) or len(sequence) != windows[i][2] - windows[i][1]:
                raise ValueError('blastdbcmd returned an unexpected window length or record count')
            out.write(f'>PTWIN{i}\n{sequence}\n')
            count += 1
    if count != len(windows):
        raise ValueError(f'blastdbcmd returned {count} of {len(windows)} requested windows')
    os.unlink(extracted)
    return subject


def local_rows(args, api, gene, windows, lengths, assembly_db, work_dir, indexed, offsets, meta, processes):
    query = os.path.join(work_dir, 'gene_queries.fa')
    with open(query, 'wb') as out:
        for qid in gene.queries:
            offset, size = offsets[qid]
            indexed.seek(offset)
            out.write(indexed.read(size))
    subject = extract_windows(args, api, assembly_db, work_dir, windows, processes)
    # One local database serves all query batches for this gene. Parallelism
    # is across genes; each gene's BLAST process uses exactly one thread.
    with tempfile.TemporaryDirectory(prefix='local_db_', dir=work_dir) as local_dir:
        local_db = os.path.join(local_dir, 'windows')
        processes.run([args.makeblastdb, '-in', subject, '-dbtype', 'nucl',
                        '-parse_seqids', '-blastdb_version', '5', '-out', local_db],
                       stdout=subprocess.DEVNULL)
        env = os.environ.copy()
        env.setdefault('BLAST_MT_QUERY_BATCH_SIZE', str(api.DEFAULT_BLAST_MT_QUERY_BATCH_SIZE))
        with closing(api.iter_exon_query_batches(query, work_dir, args.blast_query_batch_bytes)) as batches:
            for number, (batch, count, bases) in enumerate(batches, 1):
                command = api.exon_query_blast_command(args, local_db, batch)
                command[command.index('-num_threads') + 1] = '1'
                command[command.index('-word_size') + 1] = str(args.local_word_size)
                command[command.index('-evalue') + 1] = str(args.local_evalue)
                # Windows are separate subject records; never truncate a gene to
                # 100 windows merely because the coarse database had 24 contigs.
                command[command.index('-max_target_seqs') + 1] = str(max(args.max_target_seqs, len(windows)))
                if args.blast_query_batch_bytes:
                    print(f'Local BLAST query batch {number}: {count} queries, {bases} bases, '
                          f'{os.path.getsize(batch)} FASTA bytes (limit {args.blast_query_batch_bytes}); {work_dir}', file=sys.stderr)
                print('Local BLAST:', ' '.join(command), file=sys.stderr)
                proc = processes.start(command, stdout=subprocess.PIPE, text=True, env=env)
                try:
                    for line in proc.stdout:
                        fields = line.rstrip('\n').split('\t')
                        wid = int(fields[1].removeprefix('lcl|').removeprefix('PTWIN'))
                        contig, start, _ = windows[wid]
                        fields[1] = contig
                        fields[8] = str(int(fields[8]) + start)
                        fields[9] = str(int(fields[9]) + start)
                        fields[14] = str(lengths[contig])
                        # Expand only this gene's aliases, including shared queries.
                        yield from api.tabular_line_to_alias_alignments('\t'.join(fields), {}, meta, gene.queries)
                    code = proc.wait()
                    if code:
                        raise subprocess.CalledProcessError(code, command)
                finally:
                    proc.stdout.close()
                    if proc.poll() is None:
                        proc.terminate()
                        proc.wait()
                    processes.release(proc)


def write_gene_result(gid, windows, args, api, genes, lengths, assembly_db,
                      work_dir, query_path, offsets, meta, processes):
    result = os.path.join(work_dir, f'gene_{gid}.rows')
    try:
        with tempfile.TemporaryDirectory(prefix=f'gene_{gid}_', dir=work_dir) as gene_dir:
            with open(query_path, 'rb') as indexed, open(result, 'wb') as out:
                for row in local_rows(args, api, genes[gid], windows, lengths, assembly_db,
                                      gene_dir, indexed, offsets, meta, processes):
                    pickle.dump(row, out)
        return result
    except BaseException:
        if os.path.exists(result):
            os.unlink(result)
        raise


def parallel_local_rows(jobs, total, args, api, genes, lengths, assembly_db,
                        work_dir, query_path, offsets, meta):
    """Bound outstanding jobs and spool results instead of buffering HSPs in RAM."""
    workers = min(args.threads, total)
    print(f'Local realignment queue: {total} genes; up to {workers} concurrent genes; '
          '1 BLAST thread per gene', file=sys.stderr)
    processes = LocalProcesses()
    pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix='local-blast')
    pending = {}
    jobs = iter(jobs)
    completed = 0

    def submit_one():
        job = next(jobs, None)
        if job is None:
            return
        gid, windows = job
        future = pool.submit(write_gene_result, gid, windows, args, api, genes, lengths,
                             assembly_db, work_dir, query_path, offsets, meta, processes)
        pending[future] = gid

    try:
        for _ in range(workers):
            submit_one()
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            # Surface errors before starting any additional jobs.
            for future in done:
                future.result()
            for future in sorted(done, key=lambda item: pending[item]):
                gid = pending.pop(future)
                result = future.result()
                submit_one()
                with open(result, 'rb') as handle:
                    while True:
                        try:
                            row = pickle.load(handle)
                        except EOFError:
                            break
                        yield row
                os.unlink(result)
                completed += 1
                print(f'Local realignment completed {completed}/{total} genes '
                      f'(reference gene index {gid + 1})', file=sys.stderr)
    finally:
        processes.cancel()
        for future in pending:
            future.cancel()
        pool.shutdown(wait=True, cancel_futures=True)


def iter_gene_window_rows(args, lines, assembly_db, work_dir, api, seq_map, meta, alias_map):
    if not alias_map:
        raise ValueError('Two-stage alignment requires the exon alias table from the database builder')
    genes, query_genes = collect_genes(alias_map)
    path = os.path.join(work_dir, 'gene_windows.sqlite')
    with closing(sqlite3.connect(path)) as db:
        db.executescript('''
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
            PRAGMA temp_store=FILE;
            PRAGMA cache_size=-32768;
            CREATE TABLE seeds (gid INTEGER, contig TEXT, start INTEGER, end INTEGER);
            CREATE TABLE hits (gid INTEGER, contig TEXT, start INTEGER, end INTEGER,
                               strand TEXT, exon TEXT, payload BLOB);
            CREATE TABLE windows (gid INTEGER, contig TEXT, start INTEGER, end INTEGER, refine INTEGER);
        ''')
        lengths, seeded = {}, set()
        known_contigs = set(seq_map.values())
        hsps = 0
        minimap = args.candidate_aligner == 'minimap2'
        for line in lines:
            if minimap:
                qid, contig, length, start, end = parse_candidate(line)
            else:
                fields = line.rstrip('\n').split('\t')
                if len(fields) != len(api.BLAST_TABULAR_FIELDS):
                    raise ValueError('Malformed first-pass BLAST record')
                qid = fields[0]
                contig = api.replace_blast_ord_id(fields[1], seq_map)
                fields[1] = contig
                length = int(fields[14])
                left, right = sorted((int(fields[8]), int(fields[9])))
                start, end = left - 1, right
            if contig not in known_contigs:
                raise ValueError(f'Candidate search returned unknown assembly sequence ID {contig!r}')
            lengths[contig] = length
            if qid not in query_genes:
                raise ValueError(f'Candidate query {qid!r} is absent from the exon alias table')
            for gid in query_genes[qid]:
                gene = genes[gid]
                seeded.add(gid)
                db.execute('INSERT INTO seeds VALUES (?, ?, ?, ?)',
                           (gid, contig, max(0, start - gene.padding), min(length, end + gene.padding)))
                if not minimap:
                    for row in api.tabular_line_to_alias_alignments('\t'.join(fields), {}, meta, gene.queries):
                        if passes_filters(row, args):
                            db.execute('INSERT INTO hits VALUES (?, ?, ?, ?, ?, ?, ?)',
                                       (gid, contig, row.query_start, row.query_end, row.strand,
                                        row.exon_id_full, pickle.dumps(row)))
            hsps += 1
            if hsps % 10000 == 0:
                db.commit()
        db.commit()
        db.executescript('''
            CREATE INDEX seed_order ON seeds (gid, contig, start, end);
            CREATE INDEX hit_location ON hits (gid, contig, start, end);
        ''')
        skipped = refined = 0
        for gid, contig in db.execute('SELECT DISTINCT gid, contig FROM seeds ORDER BY gid, contig'):
            intervals = db.execute('SELECT start, end FROM seeds WHERE gid=? AND contig=? ORDER BY start, end', (gid, contig))
            for start, end in merged_intervals(intervals):
                refine = True
                if not minimap:
                    hits = db.execute('''SELECT exon, strand, start, end FROM hits
                        WHERE gid=? AND contig=? AND start>=? AND end<=? ORDER BY exon, strand, start, end''',
                        (gid, contig, start, end))
                    refine = not balanced_hits(genes[gid].exons, hits)
                db.execute('INSERT INTO windows VALUES (?, ?, ?, ?, ?)', (gid, contig, start, end, int(refine)))
                refined += refine
                skipped += not refine
        db.commit()
        db.execute('CREATE INDEX window_gene ON windows (gid, refine)')
        if minimap:
            print(f'First pass (minimap2): {hsps} candidate chains; {len(seeded)}/{len(genes)} gene loci seeded. '
                  f'Windows: {refined} for local BLAST; balanced-window skipping disabled for approximate candidates. '
                  f'{len(genes) - len(seeded)} gene loci have no seed; unseeded loci cannot be recovered locally.', file=sys.stderr)
        else:
            print(f'First pass (BLAST): {hsps} HSPs; {len(seeded)}/{len(genes)} gene loci seeded. '
                  f'Windows: {skipped} balanced (skip), {refined} require local realignment. '
                  f'{len(genes) - len(seeded)} gene loci have no seed; unseeded loci cannot be recovered locally.', file=sys.stderr)
        for gid, contig, start, end in db.execute('SELECT gid, contig, start, end FROM windows WHERE refine=0 ORDER BY gid, contig, start'):
            for (payload,) in db.execute('SELECT payload FROM hits WHERE gid=? AND contig=? AND start>=? AND end<=?', (gid, contig, start, end)):
                yield pickle.loads(payload)
        if not refined:
            return
        query_path = os.path.join(work_dir, 'indexed_queries.fa')
        offsets = index_queries(api, args.exon_fasta or f'{args.db}.exons.fa', query_path)
        total, = db.execute('SELECT COUNT(DISTINCT gid) FROM windows WHERE refine=1').fetchone()
        def jobs():
            for (gid,) in db.execute('SELECT DISTINCT gid FROM windows WHERE refine=1 ORDER BY gid'):
                windows = list(db.execute('SELECT contig, start, end FROM windows WHERE gid=? AND refine=1 ORDER BY contig, start', (gid,)))
                yield gid, windows
        yield from parallel_local_rows(jobs(), total, args, api, genes, lengths,
                                       assembly_db, work_dir, query_path, offsets, meta)
