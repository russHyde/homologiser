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

#' map_to_ensembl_homologues_with_biomart
#'
#' A function to map from a set of `gene_ids` in one species (sp1) to the
#'   ensembl-gene id of their homologues in a separate species (sp2). The
#'   mappings are performed using biomart databases.
#'
#'
#' @param        gene_ids      A vector of gene_ids, of format matching
#'   `idtype_sp1`.
#' @param        dataset_sp1   A `biomaRt::Mart` object for the first species.
#' @param        sp1           Ensembl contraction of the first species' name
#'   (eg, hsapiens).
#' @param        sp2           Ensembl contraction of the second species' name
#'   (eg, mmusculus).
#' @param        idtype_sp1    The type of id provided for the input species
#'   (sp1). This should be either "entrezgene" or "ensembl_gene_id".
#' @param        host          The URL for the host-site for the biomart
#'   database.
#' @param        mart_name     The name of the biomart database to be used.
#'
#' @return       A data.frame with columns ID.sp1 and
#'   ENSEMBLGENE.sp2.
#'
#' @importFrom   magrittr      %>%
#' @importFrom   dplyr         arrange_
#'
map_to_ensembl_homologues_with_biomart <- function(gene_ids = character(0),
                                                   dataset_sp1 = NULL,
                                                   sp1 = "hsapiens",
                                                   sp2 = "mmusculus",
                                                   idtype_sp1 =
                                                     "ensembl_gene_id",
                                                   host = "www.ensembl.org",
                                                   mart_name =
                                                     "ENSEMBL_MART_ENSEMBL") {
  results_empty <- data.frame(
    ID.sp1 = character(0),
    ENSEMBLGENE.sp2 = character(0),
    stringsAsFactors = FALSE
  )

  # Check that the idtype is a valid input to the function
  # - it must be either "entrezgene" or "ensembl_gene_id"
  stopifnot(idtype_sp1 %in% c("entrezgene", "ensembl_gene_id"))


  # Check that gene_ids is a vector (it can't be missing)
  # - If gene_ids is empty, return empty results
  # - If gene_ids is NULL or is not a vector, throw an error
  stopifnot(!is.null(gene_ids) && is.vector(gene_ids))
  if (length(gene_ids) == 0) {
    return(results_empty)
  }

  # The homologue column for sp2 in the biomaRt object for sp1 is obtained:
  homologue_field <- get_ensembl_homologue_field(sp2)

  # Checks on dataset_sp1
  # - If it is defined, it must
  #   a) be a mart; b) have a homologue field for sp2; and c) have a field for
  #   the input idtype_sp1
  # - If it isn't defined, use the default mart for species 1 (and check it)

  # TODO: remove this:
  # - since this function is not exported, we can ensure all callers
  #   pass a valid Mart and so don't need default-mart behaviour
  #
  # - then move mart-validity checks into select_and_filter

  dataset_sp1 <- if (is.null(dataset_sp1) && !is.null(sp1)) {
    use_default_mart(sp1, host, mart_name)
  } else {
    dataset_sp1
  }

  if (!is.null(dataset_sp1)) {
    stopifnot(
      is_valid_mart(dataset_sp1, idtype_sp1) &&
        is_valid_mart(dataset_sp1, homologue_field)
    )
  } else {
    stop("User should provide a biomart object or name for species1")
  }

  # If the input ids are non-ensembl, they have to be mapped to ensembl-ids
  #   before use in homology searches
  ensembl_ids_sp1 <- if (idtype_sp1 == "ensembl_gene_id") {
    gene_ids
  } else {
    # Map non-ensembl ids for sp1 to ensembl_gene_id
    id_to_ensembl_sp1 <- select_and_filter(
      filters = idtype_sp1, # should be entrezgene as of 2017/03/03
      values = gene_ids,
      attributes = c(idtype_sp1, "ensembl_gene_id"),
      mart = dataset_sp1
    )
    id_to_ensembl_sp1[, "ensembl_gene_id"]
  }

  # Map ensembl_gene_id for sp1 to ensembl_gene_id for sp2
  ensemblmap_sp1_sp2 <- select_and_filter(
    filters = "ensembl_gene_id",
    values = ensembl_ids_sp1,
    attributes = c("ensembl_gene_id", homologue_field),
    mart = dataset_sp1
  )

  # Map input_id for sp1 to ensembl_gene_id for sp2
  id_to_ensembl_sp1_sp2 <- if (idtype_sp1 == "ensembl_gene_id") {
    ensemblmap_sp1_sp2
  } else {
    merge(
      x = id_to_ensembl_sp1,
      y = ensemblmap_sp1_sp2,
      by = "ensembl_gene_id"
    )
  }
  # An entry should be present in the output for every gene that was in the
  #   input; even if if maps to NA
  missing_inputs <- setdiff(
    gene_ids,
    id_to_ensembl_sp1_sp2[, idtype_sp1]
  )

  # Output should be sorted by input id and then by output id, and should have
  #   column names ID.sp1 and ENSEMBLGENE.sp2
  data.frame(
    ID.sp1 = c(
      id_to_ensembl_sp1_sp2[, idtype_sp1],
      missing_inputs
    ),
    ENSEMBLGENE.sp2 = c(
      id_to_ensembl_sp1_sp2[, homologue_field],
      rep(NA, length(missing_inputs))
    ),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange_(.dots = c("ID.sp1", "ENSEMBLGENE.sp2"))
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

#' select_and_filter : extracts data from biomaRt, removes missing values
#'
#' @param        filters       Vector of filters for use in a biomaRt select().
#' @param        values        List (if multiple filters used) or vector of
#'   values for use in biomaRt::select().
#' @param        attributes    Vector of field-names for return from a biomaRt
#'   query.
#' @param        mart          A biomaRt `mart` object.
#'
#' @return       A data.frame containing columns for each of the `attributes`
#'   that were selected; and where only complete observations were kept (ie,
#'   there is no missing values in the returned data).
#'
#' @importFrom   biomaRt       getBM
#' @importFrom   stats         na.omit   setNames
#' @importFrom   magrittr      %>%
#'

select_and_filter <- function(filters,
                              values,
                              attributes,
                              mart) {
  # TODO: test this with no input `values` - it throws a warning at present

  # If none of the provided values contain valid (non-NA) values, return an
  #   empty dataframe that has the same column names as the requested
  #   attributes
  # Otherwise, search biomaRt using the values, and filter to keep only the
  #   results that are 'complete' (ie, have no NAs)

  null_result <- attributes %>%
    lapply(
      function(x) character(0)
    ) %>%
    data.frame(
      stringsAsFactors = FALSE
    ) %>%
    setNames(attributes)

  values_filtered <- na.omit(values)

  if (length(values_filtered) == 0) {
    return(null_result)
  }

  bm_result <- biomaRt::getBM(
    filters = filters,
    values = values_filtered,
    attributes = attributes,
    mart = mart
  ) %>%
    lapply(
      function(x) {
        x <- as.character(x)
        x[x == ""] <- NA
        x
      }
    ) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    na.omit()

  bm_result
}

###############################################################################
