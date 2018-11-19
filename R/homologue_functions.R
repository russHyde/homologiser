#' keep_only_one_to_one_homologues
#'
#' Takes two homologue maps, as returned by `map_to_homologues`, and drops any
#'   species1 to species2 mappings in the `forward_map` that aren't 1:1
#'   homology mappings. This occurs in two cases: i) if ID1 in species1 maps
#'   to multiple IDs in species2; or ii) for each ID1 in species1, if any of
#'   the IDs to which that gene maps (in species2) map back to more than one
#'   gene in species1.
#'
#' Starting with `forward_map`, if a non-1:1 mapping is found, no filtering
#'   takes place, and ID1 is unaffected but ID2 is set to NA.
#'
#' @param        forward_map   A `data.frame` containing colnames `ID.sp1` and
#'   `ID.sp2`.
#' @param        reverse_map   A `data.frame` containing colnames `ID.sp1` and
#'   `ID.sp2`.
#'
#' @importFrom   dplyr         filter_   select_
#'
keep_only_one_to_one_homologues <- function(
                                            forward_map,
                                            reverse_map) {
  check_validity <- function(dframe, dframe_name) {
    if (!is.data.frame(dframe) ||
      !all(c("ID.sp1", "ID.sp2") %in% colnames(dframe))
    ) {
      stop(
        paste(
          dframe_name,
          "should be a dataframe and have columns `ID.sp1` and `ID.sp2`"
        )
      )
    }
  }
  if (missing(forward_map) || missing(reverse_map)) {
    stop("dataframes forward_map & reverse_map should be defined")
  }
  check_validity(forward_map, "forward_map")
  check_validity(reverse_map, "reverse_map")

  one_to_many_a <- which_mappings_are_one_to_many(
    forward_map, "ID.sp1", "ID.sp2"
  )
  one_to_many_b <- which_mappings_are_one_to_many(
    reverse_map, "ID.sp1", "ID.sp2"
  )

  forward_map %>%
    dplyr::filter_(~!(ID.sp1 %in% one_to_many_a)) %>%
    dplyr::filter_(~!(ID.sp2 %in% one_to_many_b)) %>%
    merge(
      forward_map[, "ID.sp1", drop = FALSE],
      all.y = TRUE,
      by = "ID.sp1"
    ) %>%
    unique()
}

###############################################################################

#' which_mappings_are_one_to_many
#'
#' Takes a dataframe, and identifies any entries in a given column
#'   (`seed_col`) that have more than one associated entry in another column of
#'   that data-frame (`target_col`).
#'
#' NA-containing (seed_col, target_col) pairs are removed prior to determining
#'   the number of mappings if `na_rm=TRUE`.
#'
#' @param        x             Some data.frame.
#' @param        seed_col      A single string, one of the column names in x.
#' @param        target_col    A single string, one of the column names in x.
#' @param        na_rm         Boolean. Should NA-containing rows be removed
#'   from the (seed_col, target_col) sub-dataframe of x before calling one-many
#'   status?
#'
#' @importFrom   magrittr      %>%   extract2
#' @importFrom   dplyr         select_
#'
#' @include      data_manipulation.R

which_mappings_are_one_to_many <- function(
                                           x,
                                           seed_col,
                                           target_col,
                                           na_rm = TRUE) {
  if (missing(x)) {
    stop("x should be a defined dataframe in which_mappings_are_one_to_many")
  }
  if (missing(seed_col) || !(seed_col %in% colnames(x))) {
    stop("seed_col should be a colname of x in which_mappings_are_one_to_many")
  }
  if (missing(target_col) ||
    !(target_col %in% colnames(x)) ||
    target_col == seed_col
  ) {
    stop(
      paste(
        "`target_col` should be in the `colnames` of `x` in",
        "`which_mappings_are_one_to_many`"
      )
    )
  }

  na_dropper <- if (na_rm) {
    drop_incomplete_cases
  } else {
    identity
  }

  duplicated <- x %>%
    dplyr::select_(.dots = c(seed_col, target_col)) %>%
    na_dropper() %>%
    unique() %>%
    magrittr::extract2(seed_col) %>%
    get_duplicates()

  return(duplicated)
}
