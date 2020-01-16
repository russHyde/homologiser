#' generic for convert_species
#'
#' Given a dataset and a data-frame that contains ID-mappings from one species
#' (source) to another (target), return a version of \code{x} where the IDs
#' have been mapped to the IDs in the target species.
#'
#' @param   x   A data.frame or data.frame extension.
#' @param   homologues   A two-column data.frame, the first column containing
#'   the IDs for the source-species, the second containing the IDs from the
#'   target-species.
#' @param   ...   Further arguments
#'
#' @return   An object of the same type as \code{x}.
#'
#' @export
convert_species <- function(x, homologues, ...) {
  UseMethod("convert_species", x)
}

#' @export
convert_species.default <- function(x, homologues, ...) {
  stop("Unimplemented: default `convert_species` method")
}

#' @param   which_column   Column-name or index: which column in \code{x}
#'   contains the IDs that should be mapped to the target species? Default is
#'   to use the first column.
#' @param   new_column_name   Optional: once the IDs have been transferred from
#'   the source to the target species, what should be the name of the column of
#'   target-species IDs in the resulting data.frame? Default is the column name
#'   for the target species in the \code{homologues} data.frame.
#' @export

# nolint start
convert_species.data.frame <- function(
                                       # nolint end
                                       x,
                                       homologues,
                                       which_column = 1,
                                       new_column_name =
                                         colnames(homologues)[2],
                                       ...) {
  .get_homologues <- function(genes, hom) {
    idx <- match(genes, hom[[1]])
    hom[[2]][idx]
  }
  .filter_dataset <- function(.x) {
    keep_rows <- which(.x[[which_column]] %in% homologues[[1]])
    .x[keep_rows, ]
  }
  .update_genes_column <- function(.x, .column_index) {
    genes <- .x[[.column_index]]
    .x[[.column_index]] <- .get_homologues(genes, homologues)
    colnames(.x)[.column_index] <- new_column_name
    .x
  }

  column_index <- if (is.numeric(which_column)) {
    which_column
  } else {
    which(colnames(x) == which_column)
  }

  x_subset <- .filter_dataset(x)

  .update_genes_column(x_subset, column_index)
}
