###############################################################################

context("Tests for homology-mapping functions")

###############################################################################
## helpers
df <- function(
               ...) {
  data.frame(..., stringsAsFactors = FALSE)
}

###############################################################################

test_that("which_mappings_are_one_to_many: invalid input", {
  expect_error(
    object = which_mappings_are_one_to_many(),
    info = "No arguments to which_mappings_are_one_to_many"
  )
  expect_error(
    object = which_mappings_are_one_to_many(df(
      x = character(0),
      y = character(0)
    )),
    info = "No column arguments to which_mappings_are_one_to_many"
  )
  expect_error(
    object = which_mappings_are_one_to_many(
      df(
        x = character(0),
        y = character(0)
      ),
      target_col = "y"
    ),
    info = "No seed_col argument to which_mappings_are_one_to_many"
  )
  expect_error(
    object = which_mappings_are_one_to_many(
      df(
        x = character(0),
        y = character(0)
      ),
      seed_col = "x"
    ),
    info = "No target_col argument to which_mappings_are_one_to_many"
  )
  expect_error(
    object = which_mappings_are_one_to_many(
      df(
        x = character(0),
        y = character(0)
      ),
      seed_col = "z", target_col = "y"
    ),
    info = paste(
      "`seed_col` is absent from `colnames` of `x` in",
      "`which_mappings_are_one_to_many`"
    )
  )
  expect_error(
    object = which_mappings_are_one_to_many(
      df(
        x = character(0),
        y = character(0)
      ),
      seed_col = "x", target_col = "z"
    ),
    info = paste(
      "`target_col` is absent from `colnames` of `x` in",
      "`which_mappings_are_one_to_many`"
    )
  )
  expect_error(
    object = which_mappings_are_one_to_many(
      df(
        x = character(0),
        y = character(0)
      ),
      seed_col = "x", target_col = "x"
    ),
    info = paste(
      "`seed_col` and `target_col` should differ in",
      "`which_mappings_are_one_to_many`"
    )
  )
})

test_that("which_mappings_are_one_to_many: valid input", {
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = letters[1:3], y = 1:3),
      seed_col = "x",
      target_col = "y"
    ),
    expected = character(0),
    info = "All input map to a distinct output: no one-many"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = c("b", "b"), y = 2:3),
      seed_col = "x",
      target_col = "y"
    ),
    expected = "b",
    info = "Single one-many"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = letters[1:3], y = rep(1, 3)),
      seed_col = "x",
      target_col = "y"
    ),
    expected = character(0),
    info = "All input map to the same single output: no one-many"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = c("b", "b"), y = c(2, 2)),
      seed_col = "x",
      target_col = "y"
    ),
    expected = character(0),
    info = "Repeated one-one mapping"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = "b", y = c(2, NA)),
      seed_col = "x",
      target_col = "y",
      na_rm = TRUE
    ),
    expected = character(0),
    info = "One-one mapping, once an NA is dropped"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = "b", y = c(2, NA)),
      seed_col = "x",
      target_col = "y",
      na_rm = FALSE
    ),
    expected = "b",
    info = "Mapping to both NA and non-NA value, with na_rm = FALSE gives
1:many mapping"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = "b", y = c(2, 2), z = c(3, 4)),
      seed_col = "x",
      target_col = "y"
    ),
    expected = character(0),
    info = "1:1 mapping between cols x and y, but z column means there's
two rows"
  )
  expect_equal(
    object = which_mappings_are_one_to_many(
      df(x = "b", y = c(2, 2), z = c(3, 4)),
      seed_col = "y",
      target_col = "z"
    ),
    expected = 2,
    info = "Test with a different column name and datatype"
  )
})

###############################################################################
test_that("keep_only_one_to_one_homologues: invalid input", {
  # Input validity:
  expect_error(
    object = keep_only_one_to_one_homologues(),
    info = "No arguments to keep_only_one_to_one_homologues"
  )

  expect_error(
    object = keep_only_one_to_one_homologues(
      forward_map = "Not a dataframe",
      reverse_map = df(ID.sp1 = 1:3, ID.sp2 = 1:3)
    ),
    info = paste(
      "First arg to `keep_only_one_to_one_homologues` should be a `data.frame`"
    )
  )

  expect_error(
    object = keep_only_one_to_one_homologues(
      forward_map = df(ID.sp1 = 1:3, ID.sp2 = 1:3),
      reverse_map = "Not a dataframe"
    ),
    info = paste(
      "Second arg to `keep_only_one_to_one_homologues` should be a",
      "`data.frame`"
    )
  )

  # ID.sp1 and ID.sp2 should be the colnames of forward_map and reverse_map
  expect_error(
    object = keep_only_one_to_one_homologues(
      forward_map = df(ID.sp1 = 1:3, ID.sp2 = 1:3),
      reverse_map = df(Wrong = 1:2, COLNAMES = 3:4)
    ),
    info = "Second arg should have ID.sp1 and ID.sp2 as colnames"
  )
  expect_error(
    object = keep_only_one_to_one_homologues(
      forward_map = df(Wrong = 1:2, COLNAMES = 3:4),
      reverse_map = df(ID.sp1 = 1:3, ID.sp2 = 1:3)
    ),
    info = "First arg should have ID.sp1 and ID.sp2 as colnames"
  )
})

test_that("keep_only_one_to_one_homologues: valid input", {
  # one:one mappings in both directions:
  a <- df(ID.sp1 = 1:3, ID.sp2 = letters[1:3])
  b <- df(ID.sp1 = letters[1:3], ID.sp2 = 1:3)
  expect_equal(
    object = keep_only_one_to_one_homologues(a, b),
    expected = a,
    info = "one output for each input in both directions: no change"
  )

  # one:one mappings in both directions with named arguments:
  a <- df(ID.sp1 = 1:3, ID.sp2 = letters[1:3])
  b <- df(ID.sp1 = letters[1:3], ID.sp2 = 1:3)
  expect_equal(
    object = keep_only_one_to_one_homologues(
      forward_map = a,
      reverse_map = b
    ),
    expected = a,
    info = "one output for each input in both directions: named arguments"
  )

  # one:many mapping
  a <- df(ID.sp1 = 1, ID.sp2 = letters[1:3])
  b <- df(ID.sp1 = letters[1:3], ID.sp2 = 1)
  expect_equal(
    object = keep_only_one_to_one_homologues(a, b),
    expected = df(ID.sp1 = 1, ID.sp2 = as.character(NA)),
    info = "one-many mapping: mappings converted to a single NA"
  )

  # many:one mapping
  a <- df(ID.sp1 = 1:3, ID.sp2 = "a")
  b <- df(ID.sp1 = "a", ID.sp2 = 1:3)
  expect_equal(
    object = keep_only_one_to_one_homologues(a, b),
    expected = df(ID.sp1 = 1:3, ID.sp2 = as.character(NA)),
    info = "many-one mapping: drop all rows"
  )
})

###############################################################################
