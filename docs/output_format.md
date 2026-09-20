# ParaTyper output format (3.9.2)

## Main and fragment tables

The two user-facing outputs, `SAMPLE.transcript_calls.tsv` and `SAMPLE.pseudofragments.tsv`, share the following **44 columns in order**. Both have a header, including when empty. Empty fields represent absent information. IDs are unversioned unless explicitly noted otherwise.

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `GENE_index` | Parent index linking a gene row to its child transcripts. Unique across both tables in one run; blank in transcript-only mode. |
| 2 | `transcript_id` | Transcript ID for a real isoform; gene ID for a full-gene parent. |
| 3 | `gene_id` | Reference gene ID. Shared-exon parents join member IDs with `&`; child transcripts retain their source gene ID. |
| 4 | `gene_name` | Reference gene symbol, merged-unit name, or identical-MANE representative name. |
| 5 | `transcript_type` | Source biotype; a synthetic parent is protein-coding if any member isoform is protein-coding. |
| 6 | `ifmane_transcript` | `1` for MANE, otherwise `0`; synthetic parents are not MANE. |
| 7 | `weighted_score` | Raw chain score after applicable completeness/coding/MANE multipliers. Coding multipliers apply to real transcripts only. MANE priority is also a separate ranking preference. |
| 8 | `raw_exon_score` | Sum of distinct matched merged-interval scores minus the gene-stage insertion penalty. |
| 9 | `total_aligned_bases` | Sum of core aligned exon bases across retained exon observations; this is not unique genomic span or read depth. |
| 10 | `mean_identity` | Core exon identity percentage weighted by aligned exon bases. |
| 11 | `mean_exon_coverage` | Arithmetic mean of retained exon coverage percentages. |
| 12 | `expected_exons` | Number of eligible expected transcript exons or synthetic reference union blocks. |
| 13 | `found_exons` | Number of distinct expected exon numbers matched; inserted observations do not add completeness. |
| 14 | `missing_exons` | Expected exon count minus matched expected count. |
| 15 | `fraction_expected_found` | Matched expected count / expected count, on a 0–1 scale. |
| 16 | `found_exon_numbers` | Reference exon labels for retained observations, ordered by target position. Can include repeated/inserted observations; do not equate its list length with `found_exons`. |
| 17 | `query_id` | Original target FASTA record identifier. |
| 18 | `query_contig` | Target contig; may be parsed from a region-suffixed record name in `header-suffix` mode. |
| 19 | `query_start` | Assignment start, 0-based. Parent bounds cover owned merged intervals; ties cover the union of alternatives. |
| 20 | `query_end` | Assignment end, exclusive. |
| 21 | `strand` | Target orientation, `+` or `-`. |
| 22 | `exon_query_coordinates` | Retained target exon alignment intervals as `start-end`. |
| 23 | `exon_reference_coordinates` | Corresponding source annotated exon intervals as `contig:start-end`, converted to 0-based, half-open coordinates. These describe the annotated exon, which can exceed its aligned portion. |
| 24 | `exon_ids` | Unversioned reference exon IDs corresponding to the exon arrays. |
| 25 | `exon_coverages` | Core reference exon coverage percentages. |
| 26 | `exon_identities` | Core exon alignment identity percentages. |
| 27 | `alignment_AS` | Per-exon alignment score from the input alignment table. In the default anchored workflow, this is the core mismatch-adjusted score, not the anchored BLAST selection score. |
| 28 | `call_status` | `complete` when no eligible expected exon is missing; otherwise `partial`. Completeness does not imply expression, intact coding function, or perfect sequence identity. |
| 29 | `tie_group_id` | Identifier for a connected component of equally ranked assignments sharing merged intervals; blank for an untied call. |
| 30 | `tie_count` | Number of alternatives represented by the ambiguity component; `1` for an untied call. Not copy number. |
| 31 | `assignment_status` | `unique` or `tied` among retained candidates. |
| 32 | `pipeline_version` | Caller version, currently `3.9.2`. |
| 33 | `merged_interval_coordinates` | Distinct scored target intervals in increasing coordinate order. |
| 34 | `merged_interval_gene_scores` | Cached gene-specific similarity scores corresponding to those intervals, before bonuses. |
| 35 | `merged_interval_score_exon_ids` | Reference exons supplying those interval scores. |
| 36 | `merged_interval_score_exon_lengths` | Lengths of the supplying reference exons. |
| 37 | `transcript_exon_length` | Full annotated spliced length for real transcripts; reference exon-union length for synthetic parents. Used in same-gene isoform tie breaking. |
| 38 | `model_type` | `full_gene` or `transcript`. |
| 39 | `inserted_exons` | Sum of unique inserted reference block counts across separate insertion runs. Transcript-stage counts are zero. |
| 40 | `insertion_penalty` | Positive amount subtracted from gene raw score: `50 × inserted_exons`. |
| 41 | `inserted_exon_numbers` | Labels of inserted observations; repeats are retained. |
| 42 | `inserted_exon_query_coordinates` | Inserted merged target intervals as `start-end`. |
| 43 | `inserted_query_blocks` | Physical number of inserted merged target intervals, including repeated exon labels. |
| 44 | `insertion_run_unique_exons` | Unique exon counts per insertion run, in gene-strand order. A matched block ends a run. |

### Parent/child and ambiguity representation

A gene parent precedes its children:

```text
GENE_index    transcript_id    gene_id    ...    model_type
GENE_000001   GA               GA         ...    full_gene
GENE_000001   TA               GA         ...    transcript
```

This is a schematic excerpt, not a complete TSV row. Child rows and parents reuse exon evidence across stages. Additional partial isoforms can describe unoccupied evidence within the same copy; they do not add gene copies.

One tied assignment occupies one row. Semicolons separate alternatives in the same order as the ID fields; commas separate observations within an alternative. For example, `gene_id=GA;GB` and `found_exon_numbers=2,7;3,8` pair GA with exons 2/7 and GB with exons 3/8. Coordinate/ID/metric arrays retain repeated equal values across alternatives to preserve pairing. Other summary fields appear once if equal, and are semicolon-separated if they differ. A parser should broadcast a shared scalar where appropriate, not assume every field has `tie_count` values.

`GA&GB` is one shared-exon reference unit, not a tie. A tie between that unit and GC is `GA&GB;GC`. An identical-MANE representative retains its original ID and appends `merged` to its gene name. The representative name indicates that exact MANE sequence identity prevents distinguishing the merged genes by that evidence alone.

### Coordinates and exon arrays

All output intervals are 0-based and half-open: `100-110` spans ten bases. Exon arrays are ordered by increasing target coordinate even on the minus strand; exon numbers can therefore run backward. Keep aligned array positions together. Reference exon intervals describe the source model, while target intervals describe retained alignments. Parent bounds and scored merged intervals can extend beyond the individual exon representatives shown.

`expected_exons` counts the eligible reference exons retained by the pipeline. A complete call is therefore complete relative to that eligible model. Missing, masked, or filtered reference exons can change the denominator.

### Fragment routing

A parent and all its children are routed to `SAMPLE.pseudofragments.tsv` when the parent is partial, matches exactly one reference union block, and has no selected complete MANE child. Every alternative of a tied parent must qualify. A complete non-MANE child does not override this rule. The main and candidate tables do not duplicate these groups, and each table may have gaps in `GENE_index`.

These are **candidate gene-like fragments, not confirmed pseudogenes**. The table may contain pseudogene fragments or other partial matches. The structural rule alone does not establish biological function, retrotransposition, or the cause of the partial match; reference gene names identify the matching templates.
