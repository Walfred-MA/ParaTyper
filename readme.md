# ParaTyper

**Transcript-aware annotation of paralogs and partial gene copies in genome assemblies.**

ParaTyper transfers exon and transcript models from an annotated reference genome to a new genome assembly. It is designed for duplicated loci where assigning one gene name to an entire region can hide differences among paralogs, exon structures, and transcript isoforms.

ParaTyper identifies DNA-supported gene copies and transcript models, including paralogs, partial exon duplications/deletions, and sequence differences associated with the selected models. It also reports small gene-like fragments that can interfere with copy-number variation (CNV) analysis. Comparing these results with independent gene-level assignments can help identify **candidate gene conversion**.

The current release is **3.9.0**. Its two output tables describe reference-derived transcript models supported by genomic DNA. They do not establish RNA expression or discover arbitrary new splice isoforms. Sequence differences are summarized through exon identity, coverage, alignment scores, and exon structure; exact nucleotide alleles require sequence-level follow-up.

## 1. What ParaTyper does

- **Paralog assignment:** compares exon evidence among related genes, separates supported copies, and preserves unresolved cross-gene ties.
- **Exon structure:** reports complete and partial models, missing reference exons, and inserted/repeated exon-block evidence. These help investigate partial duplications and deletions.
- **Transcript isoforms:** selects annotated isoforms inside each assigned gene copy, with MANE priority by default and explicit links between gene and transcript rows.
- **Sequence differences:** reports exon identity, coverage, alignment scores, and corresponding reference/target coordinates for investigating mutations relative to the source reference.
- **Gene-like fragments:** separates a defined class of partial gene matches from the main call table, helping avoid counting small fragments as intact copies. These candidates may include pseudogene fragments, but are not classified as confirmed pseudogenes.
- **Gene conversion candidates:** supports comparison between the identity suggested by exon/transcript sequence and the identity suggested by whole-gene sequence or genomic context.

Missing or weak exon evidence can also reflect assembly gaps, divergent sequence, repeat masking, or alignment thresholds. A partial call is evidence to investigate, not a confirmed deletion. Copy number requires interpreting distinct loci and exon structures rather than counting output rows.

## 2. Install

Use Linux or macOS with **Python 3.9 or newer** and NCBI **BLAST+** (`blastn` and `makeblastdb`). The core Python scripts use only the standard library. On Windows, use a Linux environment such as WSL.

Clone this repository, enter its directory, then create the environment:

```bash
git clone https://github.com/Walfred-MA/ParaTyper.git
cd ParaTyper
conda env create -f environment.yml
conda activate ParaTyper
python scripts/install.py
python scripts/annotate_assemblies.py --version
```

`mamba env create -f environment.yml` is an alternative. If Python and BLAST+ are already installed, no package installation is required: run the scripts directly. `scripts/install.py` checks the scripts and executables; if BLAST+ is missing, it installs it with mamba into the active conda environment. It does not create an environment.

Check the bundled regression tests without downloading genomes or installing BLAST+:

```bash
python -m unittest discover -s test/regression -v
```

The repository contains:

```text
ParaTyper/
├── scripts/                 Pipeline and optional Snakefile
├── test/
│   ├── smn.gff3             SMN quick-start annotation
│   ├── AMY.gff3, C4.gff3, CYP2D.gff3, HPR.gff3
│   ├── LPA.gff3, MUC.gff3, NBPF.gff3, test.gff3
│   └── regression/          Tests and frozen, small alignment fixtures
├── docs/output_format.md    Complete output field reference
├── environment.yml
├── .github/ISSUE_TEMPLATE/  Gene exception report form
└── readme.md
```

The complete genome-wide GFF3 and genome FASTAs are downloaded separately. See [test/readme.md](test/readme.md) for annotation contents and provenance.

## 3. Get started: SMN on CHM13

Run the following commands from the repository root. Here **GRCh38/hg38 is the annotated source reference** and **CHM13 is the assembly being annotated**. The bundled `test/smn.gff3` uses GRCh38 coordinates; it must not be paired with CHM13 as `--reference`.

### Download the genomes and GENCODE v50 annotation

Download the **full GRCh38.p14 FASTA, including alternate loci, haplotypes, patches, and scaffolds**. Some bundled gene annotations are on these sequences. A primary-assembly-only FASTA omits relevant models. The GENCODE ALL-regions FASTA uses the same sequence names as its GFF3. [GENCODE release 50 downloads](https://www.gencodegenes.org/human/release_50.html).

We recommend downloading the matching **GENCODE v50 comprehensive annotation**, `gencode.v50.chr_patch_hapl_scaff.annotation.gff3`, which includes chromosomes, patches, haplotypes, and scaffolds. Use this annotation for whole-genome runs. The small SMN example below uses the bundled `test/smn.gff3`.

Use the NCBI T2T-CHM13v2.0 assembly for the target. This provides the `NC_0609xx.1` contig names used in the saved regression examples. [NCBI assembly GCF_009914755.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/).

```bash
mkdir -p data
curl -fL --retry 3 \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.p14.genome.fa.gz \
  -o data/GRCh38.p14.genome.fa.gz
gzip -d data/GRCh38.p14.genome.fa.gz

curl -fL --retry 3 \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz \
  -o data/gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz
gzip -d data/gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz

curl -fL --retry 3 \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
  -o data/CHM13v2.0.fa.gz
gzip -d data/CHM13v2.0.fa.gz
```

These are whole genomes. Allow space for the uncompressed FASTAs and temporary BLAST databases. Keep the reference's original sequence names and masking; the builder uses uppercase-base counts when filtering exons.

If the builder warns that a GFF3 contig is absent from the reference FASTA, those models cannot be fully evaluated. Check the FASTA release and contig names, then rebuild with the complete matching reference.

### Run ParaTyper

Create a two-column query list with one sample name and assembly path per line:

```bash
printf 'CHM13 data/CHM13v2.0.fa\n' > query_paths.txt

python scripts/annotate_assemblies.py \
  --reference data/GRCh38.p14.genome.fa \
  --gff3 test/smn.gff3 \
  --query-list query_paths.txt \
  --exon-database-dir work/smn_database \
  --output results/smn_CHM13 \
  --jobs 1 \
  --blast-threads 8 \
  --caller-threads 8
```

Paths in a query list are relative to **the query-list file**, not the shell's working directory. Quotes are supported for paths containing spaces. One multi-contig assembly belongs on one line, so all CHM13 chromosomes contribute to the same sample.

ParaTyper provides two result tables per sample:

```text
results/smn_CHM13/CHM13/CHM13.transcript_calls.tsv
results/smn_CHM13/CHM13/CHM13.pseudofragments.tsv
```

The saved SMN regression includes complete nine-exon MANE assignments for SMN2 near `NC_060929.1:70809743-70837675` on the minus strand and SMN1 near `NC_060929.1:71381874-71409804` on the plus strand (0-based, half-open). These are fixture checks, not a complete expected-results list for the broader `smn.gff3`; changing reference sequence, annotation scope, or filtering can change results.

### Other gene sets and whole-genome annotation

Replace `test/smn.gff3` with another bundled GFF3 and use separate database/output directories. The small files retain the original annotation subsets, which may include similarly named genes and alternative loci.

To annotate all GENCODE genes on CHM13, use the recommended GENCODE v50 annotation downloaded above as `--gff3`, with the same full GRCh38 FASTA and CHM13 query list:

```bash
python scripts/annotate_assemblies.py \
  --reference data/GRCh38.p14.genome.fa \
  --gff3 data/gencode.v50.chr_patch_hapl_scaff.annotation.gff3 \
  --query-list query_paths.txt \
  --exon-database-dir work/gencode_v50_database \
  --output results/gencode_v50_CHM13 \
  --blast-threads 8 --caller-threads 8
```

Both uncompressed `.gff3` and compressed `.gff3.gz` annotations are accepted. Whole-genome annotation is substantially larger than the SMN example; the bundled regression suite does not benchmark its runtime or memory requirements.

The BLAST step searches exon queries in sequential batches of at most **110 MB**
(110,000,000 FASTA bytes, including headers). Each batch uses the requested
`--blast-threads` count and the same complete assembly database. All batches
contribute to one alignment table before transcript calling. Original exon
mappings and case masking are preserved. Temporary FASTA headers contain
only exon IDs; descriptions and transcript metadata remain in the reference files.

Use `--blast-query-batch-bytes 50000000` for smaller 50 MB batches, or
`--blast-query-batch-bytes 0` to compare with the original unbatched search.
Smaller batches may take longer because they rescan the assembly. This controls
query-file size, not a hard RAM limit. A single exon record exceeding the limit
causes an error rather than being split or dropped. Batch sizes and thread counts
are logged; completed batch files are reused and removed with the temporary database.
These options affect exon-query searches, not parsing of existing alignment files.

## 4. Method and output formats

### How the pipeline works

1. **Prepare reference exons.** Read GENCODE-style GFF3 gene/transcript/exon relationships, extract exon DNA with 60 bp flanking anchors, and retain anchored sequences with at least 50 uppercase A/C/G/T bases by default. The reference FASTA must contain the GFF3's contigs. Keep exon IDs, exon numbers, gene/transcript IDs, biotypes, and MANE tags when preparing custom subsets.
2. **Define gene units.** Genes sharing more than 100 bp of exonic reference coordinates on the same contig and strand are merged transitively into shared-exon units. Separately, genes with exactly identical sets of spliced MANE DNA sequences are represented by the first gene in GFF3 order, named `<gene>merged`. This identity includes UTRs; it is not proof that all non-MANE isoforms are identical. Shared-exon groups retain their member isoforms and do not use this representative-only shortcut.
3. **Merge queries, then align to the target.** Before BLAST, overlapping core exons within the same reference gene unit, contig, and strand are replaced by one union query with flanking anchors. Partial overlaps merge transitively; flank overlap alone does not join separate exons. Every original exon ID, boundary, and transcript association is retained in the alias metadata. Each gapped BLAST hit is projected back to each original exon and its own anchors. Defaults then require anchored coverage ≥90%, identity >95%, and alignment score >50 **per original exon**. Core-exon coordinates, CIGARs, identities, and scores are recalculated from the projected alignment.
4. **Assign gene copies.** Merge overlapping target exon hits, score each interval for each gene using its longest eligible aligned reference exon, and chain synthetic full-gene exon-union models. Gene-level competition has no protein-coding preference, allowing pseudogenes to compete on the same score scale.
5. **Assign transcripts within each copy.** Real isoforms compete within the selected parent's owned intervals and reference locus. MANE has priority by default; protein-coding transcripts receive a default 10-fold score multiplier. Same-gene ties prefer longer annotated spliced transcripts, then stable IDs. Cross-gene ties remain explicit.
6. **Separate candidate fragments.** Write the main calls and the structurally defined candidate fragment table. Both files are written even when empty.

The exon similarity score is `100 × (L − 4 × (L − identical_bases)) / L`, where `L` is reference exon length. Distinct merged intervals contribute once to a chain. Gene-stage inserted runs cost 50 per unique reference exon-block number per run; a run with more than 20 unique blocks is disallowed. Skipped reference blocks have zero cost. Complete models receive a default twofold multiplier. Scores rank candidates and are not probabilities.

A `full_gene` model is the **union of annotated exons**, not an alignment of the complete genomic gene with introns. Overlapping alternative exons merge into blocks. A real transcript can therefore have more exons than its parent has union blocks. A complete synthetic parent does not mean all its alternative exons occur together in one RNA.

The pre-BLAST query merge reduces repeated searches; the later gene-level union
prevents double-counting during scoring. For a merged query, the raw BLAST command
uses `-perc_identity 0` and omits a whole-query coverage cutoff, because an intact
short exon can occupy a small part of a longer union. The configured identity and
coverage thresholds are applied to each projected original exon instead. Merged
queries can change BLAST alignment context and search statistics; exact equivalence
across all loci is not assumed.

Rebuild an existing exon database to generate the merged queries. The Python
runner accepts `--force-rebuild-database`, which also refreshes the sample results.
For standalone or Snakemake runs, rerun the database-building step before alignment.

### Output files

The user-facing outputs are **`SAMPLE.transcript_calls.tsv`** and **`SAMPLE.pseudofragments.tsv`**. Both use the same schema and are written with headers even when no calls qualify. Intermediate files and caches are managed internally by the pipeline.

| File | Contents |
| --- | --- |
| `SAMPLE.transcript_calls.tsv` | Selected gene parents and their transcript models, including complete/partial status and ambiguity. |
| `SAMPLE.pseudofragments.tsv` | Candidate gene-like fragments and their child transcripts; same columns as the main table. |

The two call tables have 44 headered, tab-separated columns. The [complete field reference](docs/output_format.md) describes every column and how to interpret the two tables. Key fields are:

| Fields | Meaning |
| --- | --- |
| `GENE_index`, `model_type` | Parent/child relationship; `full_gene` or `transcript`. |
| `transcript_id`, `gene_id`, `gene_name` | Unversioned source IDs and gene name. A parent uses the gene ID in both ID columns. |
| `query_contig`, `query_start`, `query_end`, `strand` | Target location and orientation. |
| `expected_exons`, `found_exons`, `missing_exons`, `call_status` | Reference-model completeness among eligible exons; `complete` or `partial`. |
| `exon_query_coordinates`, `exon_reference_coordinates`, `exon_ids` | Corresponding target and reference exon evidence. |
| `exon_coverages`, `exon_identities`, `alignment_AS` | Per-exon alignment metrics. |
| `weighted_score`, `raw_exon_score`, `merged_interval_*` | Chain score and the interval evidence supplying it. |
| `inserted_exons`, `inserted_query_blocks`, `insertion_*`, `inserted_exon_*` | Repeated/inserted block evidence and penalties, primarily for gene-stage interpretation. |
| `assignment_status`, `tie_count`, `tie_group_id` | Unique assignment or unresolved cross-gene alternatives. |

**Coordinates and grouping:** GFF3 input is 1-based, inclusive; output intervals are **0-based, half-open**. Within each alternative, exon arrays follow increasing target coordinate order, including on the minus strand. Commas separate exon entries; semicolons separate tied gene alternatives. An ampersand joins members of one shared-exon gene unit (`GA&GB`), while `GA;GB` denotes competing alternatives. Parse by column names and preserve these relationships.

Each parent appears before its children, sharing `GENE_index`. The index is unique across both call tables within a run, but is not a persistent identifier across different runs. Transcript rows can include residual fragments of the same parent. Neither row counts nor `tie_count` represent gene copy number. Start CNV interpretation with distinct parent loci, then inspect completeness, fragment evidence, and ties. Reference alternate loci must not be treated as additional copies in a single haploid genome.

**Candidate gene-like fragments:** a parent moves to `SAMPLE.pseudofragments.tsv` only if it is partial, has exactly one matched reference union block, and has no selected complete MANE child. For a tied parent, every alternative must qualify. Parent and children move together. The filename is a screening label: entries may include pseudogene fragments or other partial gene matches, and are not necessarily pseudogenes. This rule does not establish biological function, retrotransposition, or the cause of the partial match. A complete single-exon gene or annotated pseudogene remains in the main table.

**Mutation interpretation:** inspect `exon_identities`, `exon_coverages`, `alignment_AS`, and the exon-structure fields in the two call tables. Use `exon_reference_coordinates` and `exon_query_coordinates` to extract and realign the corresponding reference/target sequence when establishing exact alleles, orientation, or a conversion tract. Neither exon insertions nor missing exons alone establish a nucleotide-level mutation or its mechanism.

### Options, reuse, and advanced workflows

- `--no-prefer-mane` disables transcript-stage MANE priority.
- `--no-full-gene-transcripts` runs transcript-only selection, leaving `GENE_index` blank and the fragment table empty.
- `--max-chains-per-transcript` defaults to 10 per transcript/query/strand. Copy-rich loci may require increasing it. `--max-target-seqs` limits BLAST target records, not a direct copy-number threshold.
- The runner defaults to `--query-coordinate-mode local`. The standalone caller defaults to `header-suffix`, interpreting a name ending `_start_end` as a sliced-region offset. Explicitly choose `local` for ordinary assembly records.
- `--jobs` controls simultaneous assemblies; `--blast-threads` and `--caller-threads` apply per assembly. Budget CPUs and memory for concurrent jobs.
- `--blast-query-batch-bytes` controls the BLAST query FASTA batch limit (default 110,000,000 bytes). The equivalent Snakemake key is `blast_query_batch_bytes`.
- Compatible exon databases and completed sample tables can be reused. When changing query contents or alignment/calling parameters, use a new output directory or `--force-samples`; sample reuse is not a full parameter-provenance check. Use a separate database directory per annotation set.
- `--query-fasta` treats each FASTA record as a separate sample. Use `--query-list` for ordinary multi-chromosome assemblies.

Run `python scripts/annotate_assemblies.py --help` for the available pipeline options.

An optional [Snakemake configuration example](scripts/config.example.json) is provided for [scripts/Snakefile](scripts/Snakefile). Copy it to `config.json` at the repository root, edit its paths, install Snakemake separately, and run `snakemake --snakefile scripts/Snakefile --cores 8`. The unified Python runner is the quick-start workflow.

## 5. Compare with gene-level mapping and pangene

Independent gene-level and transcript-level analyses provide different evidence about a duplicated locus. [Pangene](https://github.com/lh3/pangene) builds gene graphs from **miniprot protein-to-genome alignments**; it provides a complementary coding-gene view, rather than whole-genomic-gene DNA liftover. For intron/flanking-sequence identity, also inspect an independent whole-gene DNA alignment or syntenic mapping. Pangene's protein input does not directly assess noncoding transcripts or every pseudogene fragment.

Use the same target assemblies and compatible reference gene IDs for both analyses:

1. Run ParaTyper and retain both call tables, including the parent/child links, fragment candidates, and per-exon coordinates and metrics.
2. Run the independent gene-level method. For pangene, prepare one consistent protein set with a protein-to-transcript/gene ID mapping and align it to each assembly with miniprot. Follow the upstream instructions for choosing that protein set.
3. Match assignments by target contig, strand, exon overlap, and neighboring genes. Compare full-gene identity, copy boundaries, and exon/transcript identity. Resolve naming differences, alternate loci, shared-exon units, and identical-MANE groups before interpreting disagreement.
4. Inspect discordant loci at sequence resolution. A locus with gene-A genomic context but gene-B-like exons/transcripts is a **gene-conversion candidate when both assignments are independently well supported and refer to the same copy**. Agreement between methods alone is not evidence of conversion, and an unresolved paralog tie is not sufficient.
5. Test for a localized switch in paralog-specific sequence across the candidate tract, with flanking support for the recipient locus. Check assembly gaps/collapse, copy-boundary errors, alternative annotations, and read support before concluding that conversion occurred. Similar sequences alone do not establish donor, recipient, or direction.

For example, after installing pangene and miniprot separately and preparing `proteins.faa`, their upstream workflow can be run as:

```bash
mkdir -p results/pangene
miniprot --outs=0.97 --no-cs -Iut8 data/GRCh38.p14.genome.fa proteins.faa \
  > results/pangene/GRCh38.paf
miniprot --outs=0.97 --no-cs -Iut8 data/CHM13v2.0.fa proteins.faa \
  > results/pangene/CHM13.paf
pangene results/pangene/GRCh38.paf results/pangene/CHM13.paf \
  > results/pangene/genes.gfa
```

Retain the miniprot alignments for coordinate-level comparison; the graph alone is not a conversion callset. ParaTyper does not currently automate this comparison. See the [pangene documentation](https://github.com/lh3/pangene) for protein selection and graph filtering, and the [miniprot documentation](https://github.com/lh3/miniprot) for alignment options.

## 6. Report gene exception cases

Please open an issue in this repository's **Issues** tab and choose **Gene exception case**. The [report template](.github/ISSUE_TEMPLATE/gene_exception.md) requests:

- Gene/transcript IDs, assembly accession, contig, strand, and coordinates with their convention.
- Expected versus observed behavior, including why the expected interpretation is supported.
- ParaTyper version, exact command, reference/annotation release, and relevant settings.
- A small reproducible GFF3/FASTA example or public download links, plus relevant rows from `SAMPLE.transcript_calls.tsv` and `SAMPLE.pseudofragments.tsv`.
- Independent evidence, such as another mapping method, genome-browser views, or sequence/read support.

Useful cases include incorrect paralog assignments, missed or fused copies, partial exon duplications/deletions, overlapping genes, pseudogene fragments, unresolved isoforms, and candidate gene conversion. Small examples that reproduce a problem can become regression fixtures after review.
