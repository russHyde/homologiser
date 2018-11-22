###############################################################################

context("Tests for db-connection / validity in homologiser")

###############################################################################

test_that("get_ensembl_homologue_field", {
  expect_equal(
    get_ensembl_homologue_field("hsapiens"),
    "hsapiens_homolog_ensembl_gene",
    info = "correct formatting of homology field for a given species"
  )

  expect_error(
    get_ensembl_homologue_field(),
    info = "homologue-field input should be defined"
  )

  expect_error(
    get_ensembl_homologue_field(NULL),
    info = "homologue-field input should be non-null"
  )

  expect_error(
    get_ensembl_homologue_field(123),
    info = "homologue-field is only defined for strings"
  )

  expect_error(
    get_ensembl_homologue_field("Homo sapiens"),
    info = "no spaces allowed in the homologue-field"
  )
})

###############################################################################

test_that("is_valid_mart", {
  testthat::skip_if_not(requireNamespace("mockery", quietly = TRUE))

  # mock object for use in tests
  my_mart <- data.frame(required_id = 1:3)
  attr(my_mart, "class") <- "Mart"

  with_mock(
    listAttributes = mockery::mock(
      df(name = "required_id", description = "description")
    ),
    expect_true(
      is_valid_mart(my_mart, "required_id"),
      "requested column is present in a Mart"
    ),
    .env = "biomaRt"
  )
})

###############################################################################
