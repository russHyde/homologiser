###############################################################################

context("Tests for data-manipulation in homologiser")

###############################################################################
test_that("get_duplicates", {
  expect_equal(
    object = get_duplicates(numeric(0)),
    expected = numeric(0),
    info = "No duplicates in an empty vector"
  )
  expect_error(
    object = get_duplicates(NULL),
    info = "NULL input to get_duplicates"
  )
  expect_error(
    object = get_duplicates(x = data.frame(a = 1, b = 2)),
    info = "Input to get_duplicates should be a vector"
  )
  expect_error(
    object = get_duplicates(),
    info = "No input to get_duplicates"
  )
  expect_equal(
    object = get_duplicates(c("a", "b", "c", "d")),
    expected = character(0),
    info = "Character vector with no duplicates"
  )
  expect_equal(
    object = get_duplicates(c(4, 1, 2, 2, 3, 4)),
    expected = c(2, 4),
    info = "Duplicate-containing numeric vector"
  )
})

test_that("drop_incomplete_cases", {
  a <- data.frame(x = 1:2, y = c(0, NA), z = c(NA, 1))

  expect_error(
    drop_incomplete_cases(),
    info = "No input to drop_incomplete_cases"
  )

  expect_error(
    drop_incomplete_cases(.x = "Not a dataframe"),
    info = "Input to drop_incomplete_cases should be a data.frame"
  )

  expect_error(
    drop_incomplete_cases(.x = a, .cols = "NOT A COLUMN"),
    info = "If specified, .cols should be a subset of the colnames of .xin
drop_incomplete_cases"
  )

  expect_equal(
    object = drop_incomplete_cases(a),
    expected = data.frame(x = numeric(0), y = numeric(0), z = numeric(0)),
    info = "Filtering on all columns of a df with NAs on all rows"
  )

  expect_equal(
    object = drop_incomplete_cases(a, "x"),
    expected = a,
    info = "Filtering on a complete column: no change"
  )

  expect_equal(
    object = drop_incomplete_cases(a, "y"),
    expected = a[1, ],
    info = "Filter to drop second row"
  )

  expect_equal(
    object = drop_incomplete_cases(a, "z"),
    expected = a[2, ],
    info = "Filter to drop first row"
  )
})
