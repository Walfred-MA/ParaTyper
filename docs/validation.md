# Repository preparation checks

The packaged code is the existing ParaTyper 3.9.0 pipeline, with no changes to the calling algorithms. Regression imports and fixture locations were adapted to `scripts/` and `test/regression/`. The user-facing output contract consists of `SAMPLE.transcript_calls.tsv` and `SAMPLE.pseudofragments.tsv`; the Snakefile's default targets now match those two tables.

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
