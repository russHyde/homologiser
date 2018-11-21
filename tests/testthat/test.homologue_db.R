###############################################################################

context("Tests for db-connection / validity in homologiser")

###############################################################################

# TODO: mock out all relevant tests

test_that(
  "use_default_mart", {

    # replace with skip_if_offline() when cran/testthat is updated to 2.0.1
    skip_on_travis()

    # I wasn't able to meaningfully mock-out these tests
    # - the whole purpose of this function is to get a Mart object from the
    # biomart server

    # Test that a valid Mart object is returned
    expect_is(
      use_default_mart(),
      "Mart",
      info = "test that a valid Mart object is returned"
    )

    # Test that a valid Mart object is returned - using named args
    expect_is(
      use_default_mart(
        sp = "hsapiens",
        host = "www.ensembl.org",
        mart_name = "ENSEMBL_MART_ENSEMBL"
      ),
      "Mart",
      info = "test that a valid mart is returned"
    )

    # species name, hostname and martname should all be valid in
    # use_default_mart
    expect_error(
      use_default_mart(
        sp = "NOT_A_SPECIES",
        host = "www.ensembl.org",
        mart_name = "ENSEMBL_MART_ENSEMBL"
      ),
      info = "non-standard species name"
    )

    expect_error(
      use_default_mart(
        sp = "hsapiens",
        host = "NOT_A_URL",
        mart_name = "ENSEMBL_MART_ENSEMBL"
      ),
      info = "non-standard host URL"
    )

    expect_error(
      use_default_mart(
        sp = "hsapiens",
        host = "www.ensembl.org",
        mart_name = "NOT_A_MARTNAME"
      ),
      info = "non-standard mart_name"
    )
  }
)

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
