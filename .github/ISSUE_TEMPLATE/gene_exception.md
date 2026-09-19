---
name: Gene exception case
about: Report a difficult gene, paralog, exon structure, pseudogene, or conversion candidate
title: "[Gene case] "
---

## Gene and locus

- Gene name(s) and gene/transcript IDs:
- Target assembly accession/version and haplotype:
- Contig, start, end, and strand:
- Coordinate convention (0-based half-open or 1-based inclusive):

## Expected and observed result

Describe the expected assignment or structure, the ParaTyper result, and the evidence supporting your expectation. For conversion candidates, explain the independently supported gene-level and transcript-level assignments.

## Reproduce

- ParaTyper version (`python scripts/annotate_assemblies.py --version`):
- Python and BLAST+ versions:
- Operating system:
- Reference FASTA source/accession, including alternate-locus content:
- Annotation source/release:
- Exact command and any changes from default settings:
- Was an existing database/output directory reused?

```bash
# Paste the command here.
```

## Relevant output

Include the header and relevant rows from:

- `SAMPLE.transcript_calls.tsv` (parent and child rows with the same `GENE_index`).
- `SAMPLE.pseudofragments.tsv` if relevant.

## Small example and supporting evidence

Provide a minimal GFF3/FASTA reproducer or public download links with exact intervals. Include any relevant logs, independent alignments, browser views, or read support. State whether the example can be included as a public regression fixture.
