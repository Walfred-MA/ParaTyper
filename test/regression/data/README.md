# SMN regression fixtures

`SMN.gff3` and `SMN_align_CHM13.tsv` preserve the historical bundled example.

`SMN_CHM13_gene_filter.gff3` contains the SMN1/SMN2 transcript and exon records
from `CNVs/Data/SMN.gff3`. `SMN_CHM13_gene_filter.alignments.tsv` contains the
168 corresponding exon-hit rows from the user's
`CNVs/Out/temp/CHM13/CHM13.exon_alignments.tsv`.
`SMN_CHM13_gene_filter.exon_ids.tsv` records the 84 eligible annotated exon IDs
from the shared database's `reference_exons.exon_aliases.tsv`.

Under v3.4.2's per-gene longest-exon score for each merged interval, the regression
requires a complete MANE SMN2 call at `NC_060929.1:70809743-70837675` (-) and
a complete MANE SMN1 call at `NC_060929.1:71381874-71409804` (+), using zero-based,
right-open coordinates. Each has nine exons and eight scored intervals.
The earlier maximum-score rule tied SMN1/SMN2 at the positive locus because
short perfect exons masked sequence differences in longer exons. These checks
validate this supplied example, not an independent biological ground truth.

Equal-ranked isoforms of the same gene now choose the longer annotated spliced
transcript; output ambiguity is restricted to different genes.

## C4 tandem-copy regression

`C4_CHM13_tandem.gff3` retains the C4A/C4B annotation rows from `Data/C4.gff3`.
`C4_CHM13_tandem.alignments.tsv` contains the 538 corresponding saved CHM13
alignment rows; `C4_CHM13_tandem.exon_ids.tsv` freezes the 269 eligible exon IDs
from the shared database. The source database was built with `Data/test.gff3`,
which includes these genes plus others; the fixture restricts the annotation
and evidence to C4A/C4B.

Under v3.5.0, a C4A chain crossed from exon 34 at 31852511–31852602 to exon 35
at 31885430–31885505 and consumed the end of the second copy. Version 3.5.1
requires two complete 41-exon assignments: C4A at 31835262–31855887 and C4B
at 31868000–31888625 on NC_060930.1 (+), with identical serial/parallel output.
Coordinates are zero-based and half-open. This is a regression against the
observed chaining failure, not an independent biological identity benchmark.

## C4orf50 exon-union regression

`C4orf50_CHM13.gff3` retains the gene/transcript/exon rows for C4orf50 from
`Data/C4.gff3`. The companion alignment and exon-ID files freeze its 39 saved
CHM13 exon hits and 39 eligible IDs. Three annotated isoforms contribute 35
reference union blocks. The default full-gene candidate must produce one
complete 35-block call at `NC_060928.1:5869135-6174552` (-); disabling synthetic
models reproduces the MANE call and two partial residual isoforms.

The C4 and SMN exact MANE-coordinate regressions explicitly use
`--no-full-gene-transcripts`; additional default-mode tests verify that
synthetic models preserve gene identity and copy separation.

## AMY pseudogene regression

`AMY_CHM13.gff3` and `AMY_CHM13.alignments.tsv` preserve the user's AMY
annotation and all 626 saved exon alignments from `Out3/temp/CHM13/`.
`AMY_CHM13.exon_ids.tsv` freezes 77 eligible exon IDs; the accompanying
`AMY_CHM13.identical_paralogs.tsv` preserves the AMY1A/B/C merge.
The AMY metadata was rebuilt separately with the original reference and
builder defaults after the shared database had been reused for SMN. It
reproduced the user's original v3.7.1 calls byte-for-byte. Source hashes are
in `AMY_CHM13.provenance.json`.

Gene-stage coding bonus previously favored a partial AMY2A chain scoring
6982.683983 over a complete AMYP1 chain scoring 1200 at
103575109–103581258 (-), despite perfect alignment of all seven AMYP1 exons.
Version 3.7.2 uses a gene-stage coding multiplier of one, preserving the
transcript-stage option. Tests require three complete AMYP1 models, the
user-specified locus coordinates, preserved coding copies, no long AMY2A
parent spanning this AMYP1 locus, and identical serial/parallel output.
They also vary the coding multiplier to verify that it changes transcript
scores without changing gene identity or gene-stage scores.

## MUC20 shared-exon gene-unit regression

`MUC20_shared_exons.gff3` retains the primary chr3 MUC20/MUC20-OT1 gene,
transcript and exon records from `Data/MUC.gff3`, with caller-relevant
attributes. It contains 258 real isoforms. The two genes share 1,851 unique
exonic reference bases and form one 48-block unit, `MUC20&MUC20-OT1`.

`MUC20_CHM13_target.alignments.tsv` and `.exon_ids.tsv` freeze a targeted
BLAST run on `NC_060927.1:198600000-198770000`, using the normal 60 bp anchors,
coverage/identity/score filters and GRCh38 reference. The FASTA query name
encodes this offset; caller default header-suffix coordinates restore the
assembly positions. Provenance JSON files document inputs and hashes.

The regression requires one 46/48-block merged parent at
198625730–198709445 (-) and a MUC20 MANE child with 3/4 exons under that same
parent. It does not assert transcript completeness: the MANE second exon is
still missing from this alignment evidence. This is a targeted implementation
check, not a whole-assembly rerun or independent biological benchmark.
