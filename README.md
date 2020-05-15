<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/russHyde/homologiser.svg?branch=master)](https://travis-ci.org/russHyde/homologiser)

[![Coverage
Status](https://img.shields.io/codecov/c/github/russHyde/homologiser/master.svg)](https://codecov.io/github/russHyde/homologiser?branch=master)

homologiser
===========

{homologiser} provides a simple route to obtain to obtain the gene IDs
for the homologues of a set of genes. It calls {biomaRt}. It provides
one exported function (`map_to_homologues`) that should be used for
cross-species ID mapping.

But using {biomaRt} is not that difficult, so why write a wrapper around
it that has a really narrow purpose?

1.  {homologiser} is simpler

2.  you might (I frequently do) need to disregard those genes that have
    no homologues, or multiple homologues in the target species, or that
    map to a gene that has multiple homologues in the source species: to
    restrict to just these genes, you can use
    `map_to_homologues(blah, blah, ..., one_to_one = TRUE)`.

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
# Let's face it, you're all pinning your analyses to a specific database
# for reproducibility, aren't you?
ensembl_v84 <- "http://Mar2016.archive.ensembl.org"
human_biomart <- biomaRt::useMart(
  biomart = "ensembl", host = ensembl_v84, dataset = "hsapiens_gene_ensembl"
)
mouse_biomart <- biomaRt::useMart(
  biomart = "ensembl", host = ensembl_v84, dataset = "mmusculus_gene_ensembl"
)
```

``` r
# A selection of human IDs
human_genes <- c("ENSG00000134294", "ENSG00000284192", "ENSG00000002726")

# Get the ensembl IDs for mouse homologues:
homologiser::map_to_homologues(
  gene_ids = human_genes,
  dataset_sp1 = human_biomart, sp1 = "hsapiens", idtype_sp1 = "ensembl_gene_id",
  dataset_sp2 = mouse_biomart, sp2 = "mmusculus", idtype_sp2 = "ensembl_gene_id"
)
#> Cache found
#>            id_sp1             id_sp2
#> 1 ENSG00000002726 ENSMUSG00000029811
#> 2 ENSG00000002726 ENSMUSG00000029813
#> 3 ENSG00000002726 ENSMUSG00000039215
#> 4 ENSG00000002726 ENSMUSG00000068536
#> 5 ENSG00000134294 ENSMUSG00000022462
#> 6 ENSG00000284192               <NA>
```

Note that - “ENSGxxxxxx2726” maps to several mouse genes -
“ENSGxxxx134294” maps to a single mouse gene - “ENSGxxxx284192” maps to
no mouse genes - although that function-call was a bit wordy, human -\>
mouse is the default direction and ensembl-gene is the default ID-type,
so we only really needed to type
`map_to_homologues(human_genes, human_biomart, mouse_biomart)`

What if we only want to consider those homologue pairings where there is
a single human gene mapping to/from a single mouse gene:

``` r
# Get the ensembl IDs for mouse homologues:
homologiser::map_to_homologues(
  gene_ids = human_genes,
  dataset_sp1 = human_biomart,
  dataset_sp2 = mouse_biomart,
  one_to_one = TRUE
)
#> Cache found
#> Cache found
#>            id_sp1             id_sp2
#> 1 ENSG00000002726               <NA>
#> 5 ENSG00000134294 ENSMUSG00000022462
#> 6 ENSG00000284192               <NA>
```

Now, since x2726 (one-to-many) and x284192 (one-to-zero) don’t map
one-to-one, they have a missing value in the returned data-frame.

What if a gene is part of a set of genes that map many-to-one? For
example, this is one of the mouse genes that “ENSGxx2726” maps to:

``` r
mouse_gene <- "ENSMUSG00000029811"

homologiser::map_to_homologues(
  gene_ids = mouse_gene,
  dataset_sp1 = mouse_biomart, sp1 = "mmusculus",
  dataset_sp2 = human_biomart, sp2 = "hsapiens",
  one_to_one = FALSE
)
#> Cache found
#>               id_sp1          id_sp2
#> 1 ENSMUSG00000029811 ENSG00000002726

homologiser::map_to_homologues(
  gene_ids = mouse_gene,
  dataset_sp1 = mouse_biomart, sp1 = "mmusculus",
  dataset_sp2 = human_biomart, sp2 = "hsapiens",
  one_to_one = TRUE
)
#> Cache found
#> Cache found
#>               id_sp1 id_sp2
#> 1 ENSMUSG00000029811   <NA>
```

Since that mouse gene is part of a many-to-one mapping, it does not have
any homology partners when we restrict to one-to-one mappings (but it’s
human homologue is included when we are less strict).

Note that you can use either “ensembl\_gene\_id” or “entrezgene” as the
“idtype”

``` r
# ensembl-mouse to entrez-human
homologiser::map_to_homologues(
  gene_ids = mouse_gene,
  dataset_sp1 = mouse_biomart, sp1 = "mmusculus",
  dataset_sp2 = human_biomart, sp2 = "hsapiens", idtype_sp2 = "entrezgene",
  one_to_one = FALSE
)
#> Cache found
#> Cache found
#>               id_sp1 id_sp2
#> 1 ENSMUSG00000029811     26

# entrez-human to ensembl-mouse
homologiser::map_to_homologues(
  gene_ids = c("10000", "1234"), # AKT3 and CCR5
  dataset_sp1 = human_biomart, idtype_sp1 = "entrezgene",
  dataset_sp2 = mouse_biomart,
  one_to_one = TRUE
)
#> Cache found
#> Cache found
#> Cache found
#> Cache found
#>   id_sp1             id_sp2
#> 1  10000 ENSMUSG00000019699
#> 2   1234               <NA>
```
