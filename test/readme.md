# Example annotations and regression fixtures

The nine top-level GFF3 files are the existing small GRCh38 annotation subsets from the development data directory. They were copied without changing their records; `SMN.gff3` was renamed to **`smn.gff3`** for the quick start. The source complete annotation identifies itself as GENCODE human v50 (Ensembl 116). The complete `gencode.v50.chr_patch_hapl_scaff.annotation.gff3` is deliberately excluded.

Use the full [GENCODE GRCh38.p14 ALL-regions genome FASTA](https://www.gencodegenes.org/human/release_50.html) with these annotations. The filenames describe the original selections; some contain genes with matching name substrings and alternative-locus records beyond the central family of interest. These are workflow examples, not independently validated copy-number truth sets.

| File | Gene records | Transcript records | Notes |
| --- | ---: | ---: | --- |
| `smn.gff3` | 10 | 109 | SMN1, SMN2, their antisense annotations, and SMNDC1; includes alternate/patch loci. |
| `AMY.gff3` | 6 | 45 | AMY1A/B/C, AMY2A/B, and AMYP1. |
| `C4.gff3` | 37 | 363 | C4A/B, associated annotations, C4orf genes, and other matching names; includes alternate loci. |
| `CYP2D.gff3` | 18 | 116 | CYP2D6, CYP2D7, CYP2D8P across primary and alternate/patch loci. |
| `HPR.gff3` | 4 | 94 | GRHPR, HPRT1P1, HPRT1P2, and SHPRH in the supplied selection. |
| `LPA.gff3` | 10 | 177 | LPA, LPAL2, LPAR1–6, and HELLPAR. |
| `MUC.gff3` | 67 | 1,627 | Mucin-related names, overlapping annotations, pseudogenes, and alternate loci. |
| `NBPF.gff3` | 28 | 240 | NBPF family and pseudogenes, including alternate/patch loci. |
| `test.gff3` | 138 | 2,328 | Existing combined subset; not the complete annotation or necessarily the union of all other files. |

Counts refer to records, so separate reference loci may share a gene name. [gff3_manifest.tsv](gff3_manifest.tsv) records source-relative paths, sizes, and SHA-256 hashes for all bundled GFF3 files, including the six regression annotations. Original data content is preserved.

## Regression suite

From the repository root:

```bash
python -m unittest discover -s test/regression -v
```

The 98 existing tests cover reference exon unions, hierarchical gene/isoform calls, copy chaining, MANE preference, shared-exon groups, identical-MANE paralogs, gene ties, fragment routing, output schemas, BED conversion, database reuse, and serial/parallel consistency. Saved alignment fixtures allow testing without BLAST+ or whole-genome FASTAs. Runner tests substitute fixture alignments for the BLAST stage and execute the real caller.

The [regression fixture notes](regression/data/README.md) describe the SMN, C4, C4orf50, AMY, and MUC20 cases. Paths in their provenance files identify the original development inputs; they are historical labels, not runtime dependencies. Absolute personal-directory prefixes have been removed from JSON provenance without changing sequence/alignment hashes.

Run the full SMN example in the [main README](../readme.md) to test reference extraction and BLAST against a downloaded CHM13 assembly. Its user-facing results are `SAMPLE.transcript_calls.tsv` and `SAMPLE.pseudofragments.tsv`.
