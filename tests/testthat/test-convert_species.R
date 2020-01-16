###############################################################################

context(
  "Test functions for converting a data-structure from one species to another"
)

###############################################################################

data1 <- .tbl(
  ensembl_gene_id = c(
    "ENSMUSG00123456789",
    "ENSMUSG00987654321",
    "ENSMUSG00011122233",
    NA,
    "ENSMUSG00123456789"
  ),
  some_data = c(
    "ABC",
    "DEF",
    NA,
    "GHI",
    "JKL"
  )
)

homologues1 <- .tbl(
  ensembl_mm = c(
    "ENSMUSG00123456789", "ENSMUSG00011122233", "ENSMUSG44455566677"
  ),
  ensembl_hs = c(
    "ENSG00100100101", "ENSG00200200202", "ENSG00300300303"
  )
)

converted <- .tbl(
  ensembl_hs = c("ENSG00100100101", "ENSG00200200202", "ENSG00100100101"),
  some_data = c("ABC", NA, "JKL")
)

###############################################################################

test_that("`data.frame`s can be converted from one species to another", {
  expect_silent(
    convert_species(
      x = .tbl(
        ensembl_gene_id = character(0), some_data = character(0)
      ),
      homologues = .tbl(
        ensembl_mm = character(0), ensembl_hs = character(0)
      ),
      which_column = "ensembl_gene_id",
      new_column_name = "ensembl_hs"
    )
  )

  expect_equal(
    convert_species(
      x = data1, homologues = homologues1, which_column = "ensembl_gene_id"
    ),
    expected = converted,
    info = "when `new_column_name` is missing, use colname of homologues[, 2]"
  )

  expect_equal(
    convert_species(
      x = data1, homologues = homologues1, new_column_name = "ensembl_hs"
    ),
    expected = converted,
    info = "when `which_column` is missing, use x[, 1]"
  )

  expect_equal(
    convert_species(
      x = data1, homologues = homologues1, which_column = "ensembl_gene_id",
      new_column_name = "my_new_column"
    ),
    expected = .tbl(
      "my_new_column" = converted[["ensembl_hs"]],
      "some_data" = converted[["some_data"]]
    ),
    info = "user can specify the col-name for the homologue column"
  )

  expect_equal(
    convert_species(
      x = data1[c("some_data", "ensembl_gene_id")],
      homologues = homologues1, which_column = "ensembl_gene_id"
    ),
    expected = converted[c("some_data", "ensembl_hs")],
    info = "gene-id column can be anywhere in the input data-frame"
  )

  expect_equal(
    convert_species(
      x = data1[c("some_data", "ensembl_gene_id")],
      homologues = homologues1, which_column = 2
    ),
    expected = converted[c("some_data", "ensembl_hs")],
    info = "gene-id column can be specified by index"
  )
})

test_that("convert_species for a data.frame: input data is valid", {
  expect_error(
    convert_species(
      x = "NOT A DATA FRAME", homologues = homologues1
    ),
    regexp = "Unimplemented: default `convert_species` method",
    info = "convert_species is currently implemented for `x: data.frame` only"
  )

  # TODO

  # homologues dataset should have two columns

  # which_column should be present in the column names (if non-numeric)

  # which_column should lie in 1..ncol(x) if numeric
})
