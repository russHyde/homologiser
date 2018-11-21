###############################################################################

#' use_default_mart
#'
#' Sets up a default biomaRt `mart` object for use in searches
#'
#' @param        sp            A species name (eg, "hsapiens", "mmusculus") as
#'   used in biomaRt.
#' @param        host          The name of the site that hosts the mart
#'   database.
#' @param        mart_name     The name of the biomart dataset that should be
#'   used.
#'
#' @importFrom   RCurl         url.exists
#' @importFrom   biomaRt       useMart   useDataset

use_default_mart <- function(sp = "hsapiens",
                             host = "www.ensembl.org",
                             mart_name = "ENSEMBL_MART_ENSEMBL") {
  # returns a default biomart dataset for the given species
  stopifnot(RCurl::url.exists(host))

  biomaRt::useMart(
    biomart = mart_name,
    dataset = paste(sp, "gene", "ensembl", sep = "_"),
    host = host
  )
}

###############################################################################

#' is_valid_mart
#'
#' Checks if a user-provided object is a valid `mart` object from the biomaRt
#'   package; and if it is, where the user-specified id-type is present as one
#'   of the attributes (ie, columns) of that mart
#'
#' @param        mart          A `mart` object from the `biomaRt` package.
#' @param        id_type       Check that this identifier type is a named
#'   field / column in the `mart` object.
#'
#' @importFrom   biomaRt       listAttributes
#' @importFrom   methods       is

is_valid_mart <- function(mart,
                          id_type = "entrezgene") {
  # checks if the input is a valid biomart dataset
  # checks that the user-specified id-type is present in the attributes

  !is.null(mart) &&
    methods::is(mart, "Mart") &&
    id_type %in% biomaRt::listAttributes(mart)[, 1]
}

###############################################################################

#' get_ensembl_homologue_field
#'
#' In an ensembl biomart object, the column that contains the ensembl-gene ids
#'   for homologues in species `sp` is <sp>_homolog_ensembl_gene, where <sp>
#'   is the interpolated value of `sp` (eg, mmusculus).
#'
#' @param        sp            A species-name contraction, as used in
#'   ensembl / `biomaRt`, for example, mmusculus for mouse and hsapiens for
#'   human.
#'
get_ensembl_homologue_field <- function(sp) {
  if (missing(sp) || is.null(sp)) {
    stop("`sp` should be defined and non-null")
  }
  if (!is.character(sp) || grepl(" ", sp)) {
    stop("`sp` should be a string containing no spaces")
  }

  paste(sp, "homolog", "ensembl", "gene", sep = "_")
}

###############################################################################
