# Repository preparation checks

Version 3.9.1 corrects full-gene chaining across disconnected target loci. All
133 regression tests passed with Python 3.10.21 in both the public repository
and the local transcript source tree. New cases cover both strands, touching
versus disconnected padded windows, insertion scoring within a locus, separate
gene parents, and reference spans containing unaligned exons.

Replaying 47,552 saved NBPF/CHM13 exon alignments reproduced the old NBPF20 call
at `NC_060925.1:16173421-144542504` (128,369,083 bp). With the correction, the
NBPF20 call near 144.5 Mb is `144414626-144542504` (127,878 bp), with 123 of 143
reference union blocks found. Distant evidence competes in separate loci.
The largest selected full-gene span in the replay is 228,660 bp. The replay
produced 123 main rows and 128 fragment rows. Serial execution using the local
source and two-process execution using the public source produced byte-identical
tables. No BLAST rerun or new biological validation of paralog assignments was
performed for this correction.

The initial repository packaging used the ParaTyper 3.9.0 pipeline without
changing its calling algorithms. Regression imports and fixture locations were
adapted to `scripts/` and `test/regression/`. The user-facing output contract
consists of `SAMPLE.transcript_calls.tsv` and `SAMPLE.pseudofragments.tsv`; the
Snakefile's default targets match those two tables.

Checks performed during repository preparation:

- All **98 existing regression tests passed** with Python 3.13.1.
- After renaming the fragment output to `SAMPLE.pseudofragments.tsv`, all 98 tests passed again, including named-sample and combined-FASTA output/resume checks. Replaying the saved SMN alignments produced both tables with byte-identical contents to the earlier results; only the fragment filename changed.
- The installer/dependency check passed with BLAST+ 2.17.0+ available.
- A fresh database build, whole-CHM13 BLAST alignment, hierarchical call, and BED conversion completed with Python 3.10.21 and BLAST+ 2.17.0+.
- That SMN run produced 18 main assignment rows and a header-only fragment table, each with the current 44-column schema. Complete unique nine-exon MANE models were recovered for SMN1 at `NC_060929.1:71381874-71409804` (+) and SMN2 at `NC_060929.1:70809743-70837675` (-), using 0-based, half-open coordinates.
- The bundled-fixture SMN HTML/TSV audit and BED conversion ran successfully.
- All 15 bundled GFF3 files match their original source bytes and the hashes in `test/gff3_manifest.tsv`.
- Local documentation links resolve. The three reference/annotation/assembly download URLs in the README returned HTTP 200 to metadata requests; those large downloads were not repeated for validation.

The fresh SMN run used the available local `GRCh38_full_analysis_set_plus_decoy_hla.fa` as source and `GCF_009914755.1_T2T-CHM13v2.0_genomic.fa` as target. **That local source FASTA lacks `GL339449.2` and `KI270897.1`**, and the builder reported both omissions. This checks execution and the retained primary-locus SMN calls; it does not validate alternate-locus results from the full GENCODE GRCh38.p14 download recommended in the README. No new whole-annotation run was performed.

Generated databases, full-genome results, and audit reports are excluded from Git. The GitHub workflow is configured to run the regression suite on Python 3.9, 3.11, and 3.13 after upload; those remote jobs have not run as part of this local preparation.
