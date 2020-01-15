###############################################################################

#' drop_incomplete_cases
#'
#' Filters an input data.frame, keeping only those rows that have complete data
#' (ie no NA values) within the specified columns
#'
#' @param        .x            A data.frame
#' @param        .cols         A vector subset of the column-names/indices of
#' \code{.x}.
#'
#' @importFrom   stats         complete.cases
#' @export

drop_incomplete_cases <- function(.x,
                                  .cols = colnames(.x)) {
  if (missing(.x) || !is.data.frame(.x)) {
    stop("Dataframe .x should be defined in drop_incomplete_cases")
  }

  check_numeric <- is.numeric(.cols) && all(.cols %in% seq_along(.x))
  check_names <- is.character(.cols) && all(.cols %in% colnames(.x))
  if (!check_numeric && !check_names) {
    stop(
      paste(
        "If defined, .cols should be a subset of either the colnames or",
        "the column-indices of .x"
      )
    )
  }

  .x[complete.cases(.x[, .cols]), ]
}

###############################################################################

#' get_duplicates
#'
#' Determines which entries in an input vector are duplicated within that
#' vector and returns those __values__.
#'
#' @param        x             A vector of some form.
#'
#' @return       A vector contain the entries of `x` that are duplicated within
#' `x`.  The order of `x` is not respected and the returned vector is
#' typically smaller than `x`.
#'
#' @export

get_duplicates <- function(x) {
  if (missing(x)) {
    stop("No input to get_duplicates")
  }
  if (!is.vector(x)) {
    stop("Input to get_duplicates should be a vector")
  }
  unique(
    x[duplicated(x)]
  )
}
