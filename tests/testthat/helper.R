###############################################################################

## helpers for test_that

df <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}

###############################################################################

# Mock biomaRt objects for use in unit-testing

mock_mart <- structure(.Data = "mock-mart", class = "Mart")

mock_attribs <- df(
  name = c(
    "ensembl_gene_id", "entrezgene", "mmusculus_homolog_ensembl_gene",
    "hsapiens_homolog_ensembl_gene"
  ),
  descriptions = c(
    "Gene ID", "External Gene ID", "Mouse homologues", "Human homologues"
  )
)

###############################################################################
