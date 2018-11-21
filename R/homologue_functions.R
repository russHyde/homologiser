#' map_to_homologues
#'
#' Takes gene ids from one species and maps them to gene ids in another
#' species. Can convert entrez ids and ensembl ids.
#'
#' Note: This function is used in project `dtg_rnaseq`. But the version there
#' uses dot-separated argument names. The argnames have been converted to
#' underscore-separated. If this function-script is copied back into the
#' `dtg_rnaseq` project you will have to double-check the function-calls.
#'
#' @param        gene_ids      A vector of gene identifiers for species `sp1`
#'   and of type `idtype_sp1`.
#' @param        dataset_sp1   The name of the `biomaRt` dataset that houses
#'   the data for species `sp1`.
#' @param        dataset_sp2   The name of the `biomaRt` dataset that houses
#'   the data for species `sp2`.
#' @param        sp1           The name of species 1 (using ensembl prefixes
#'   like `hsapiens` / `mmusculus`).
#' @param        sp2           The name of species 2 (using ensembl prefixes
#'   like `hsapiens` / `mmusculus`).
#' @param        idtype_sp1    The type of identifier used for species `sp1`;
#'   this should match one of the field / column names in the biomaRt dataset
#'   for species `sp1`.
#' @param        idtype_sp2    The type of identifier used for species `sp2`;
#'   this should match one of the field / column names in the biomaRt dataset
#'   for species `sp2`.
#' @param        host          The URL for the website that hosts the biomaRt
#'   dataset.
#' @param        mart_name     The name of the biomaRt database that is to be
#'   used.
#' @param        one_to_one    Boolean. Should the function only return
#'   one-to-one homology mappings?
#'
#' @return       A data.frame with columns ID.sp1 and ID.sp2.
#'
#' @export
#'

map_to_homologues <- function(gene_ids = character(0),
                              dataset_sp1 = NULL,
                              dataset_sp2 = NULL,
                              sp1 = "hsapiens",
                              sp2 = "mmusculus",
                              idtype_sp1 = "ensembl_gene_id",
                              idtype_sp2 = "ensembl_gene_id",
                              host = "www.ensembl.org",
                              mart_name = "ENSEMBL_MART_ENSEMBL",
                              one_to_one = FALSE) {
  # Validity checks
  # If the user hasn't provided a biomaRt 'Mart' object for either species 1 or
  #   species 2, use the "<species>_gene_ensembl" mart from ensembl as default
  #
  # Then check that the provided / constructed biomart datasets are valid
  #
  if (is.null(dataset_sp1) && !is.null(sp1)) {
    dataset_sp1 <- use_default_mart(sp1)
  }
  if (is.null(dataset_sp2) && !is.null(sp2)) {
    dataset_sp2 <- use_default_mart(sp2)
  }
  stopifnot(is_valid_mart(dataset_sp1, idtype_sp1))
  stopifnot(is_valid_mart(dataset_sp2, idtype_sp2))

  stopifnot(idtype_sp1 %in% c("ensembl_gene_id", "entrezgene"))
  stopifnot(idtype_sp2 %in% c("ensembl_gene_id", "entrezgene"))

  # Check the validity of the input gene_ids and return an empty result if
  #   no gene_ids were provided
  #
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    null_results <- data.frame(
      ID.sp1 = character(0),
      ID.sp2 = character(0),
      stringsAsFactors = FALSE
    )
    return(null_results)
  }
  stopifnot(is.character(gene_ids))

  # Check that the one_to_one argument is TRUE/FALSE
  stopifnot(is.logical(one_to_one) && length(one_to_one) == 1)

  # Map from ids in species 1 to ids in species 2
  forward_map <- map_to_homologues_oneway(
    gene_ids = gene_ids,
    dataset_sp1 = dataset_sp1,
    dataset_sp2 = dataset_sp2,
    sp1 = sp1,
    sp2 = sp2,
    idtype_sp1 = idtype_sp1,
    idtype_sp2 = idtype_sp2
  )

  # If the user wants only homologues that are one:one from one species to the
  # other, we have to compare the mappings from sp1->sp2 and the mappings from
  # sp2->sp1, keeping only the genes in sp1 that map to a single gene in sp2
  # and the genes in sp2 that map to a single gene in sp1
  if (one_to_one) {
    # Make mapping from second species back to first and
    # Filter the forward map based on 1:many maps in both forward and reverse
    # map
    reverse_map <- map_to_homologues_oneway(
      gene_ids = forward_map %>%
        drop_incomplete_cases() %>%
        magrittr::extract2("id_sp2") %>%
        unique() %>%
        sort(),
      dataset_sp1 = dataset_sp2,
      dataset_sp2 = dataset_sp1,
      sp1 = sp2,
      sp2 = sp1,
      idtype_sp1 = idtype_sp2,
      idtype_sp2 = idtype_sp1
    )
    return(
      keep_only_one_to_one_homologues(forward_map, reverse_map)
    )
  } else {
    return(
      forward_map
    )
  }
}

###############################################################################

#' map_to_homologues_oneway
#'
#' Takes gene ids from one species and maps them to gene ids in another
#' species. Can convert entrez ids and ensembl ids
#'
#'  --- Non-exported function ---
#'  --- All inputs should have been validity-checked by map_to_homologues ---
#'
#'  Selects all genes from the biomart dataset for species1 that have
#'    homologues in species 2
#'  1) Map input-ids in species 1 to ensembl gene in species 2
#'  2) Map ensembl gene in species 2 to output-ids in species 2
#'
#'  Map species 1 gene_ids to ensembl-ids in species 2
#'  - Since dataset_sp1 is the mart for species 1, we don't have to pass the
#'  sp1, host or mart.name arguments into map_to_ensembl_homologues...
#'
#' @param        gene_ids      A vector of gene identifiers for species `sp1`
#'   and of type `idtype_sp1`.
#' @param        dataset_sp1   The `biomaRt` dataset that houses the data for
#'   species `sp1`.
#' @param        dataset_sp2   The `biomaRt` dataset that houses the data for
#'   species `sp2`.
#' @param        sp1           The name of species 1 (using ensembl prefixes
#'   like `hsapiens` / `mmusculus`).
#' @param        sp2           The name of species 2 (using ensembl prefixes
#'   like `hsapiens` / `mmusculus`).
#' @param        idtype_sp1    The type of identifier used for species `sp1`;
#'   this should match one of the field / column names in the biomaRt dataset
#'   for species `sp1`.
#' @param        idtype_sp2    The type of identifier used for species `sp2`;
#'   this should match one of the field / column names in the biomaRt dataset
#'   for species `sp2`.
#'
#' @importFrom   dplyr         arrange_   bind_rows   filter_   select_
#' @importFrom   magrittr      %>%
#'

map_to_homologues_oneway <- function(gene_ids = character(0),
                                     dataset_sp1 = NULL,
                                     dataset_sp2 = NULL,
                                     sp1 = "hsapiens",
                                     sp2 = "mmusculus",
                                     idtype_sp1 = "ensembl_gene_id",
                                     idtype_sp2 = "ensembl_gene_id") {
  sp1_to_sp2_ensembl <- map_to_ensembl_homologues_with_biomart(
    gene_ids = gene_ids,
    dataset_sp1 = dataset_sp1,
    sp2 = sp2,
    idtype_sp1 = idtype_sp1
  )

  # Map species 2 ensembl-ids to species 2 output idtype
  sp2_ensembl_to_output <- if (idtype_sp2 == "ensembl_gene_id") {
    ens_sp2 <- unique(sp1_to_sp2_ensembl[, "ensembl_gene_sp2"])
    data.frame(
      ensembl_gene_id = ens_sp2,
      id_sp2 = ens_sp2,
      stringsAsFactors = FALSE
    )
  } else {
    select_and_filter(
      filters = "ensembl_gene_id",
      values = sp1_to_sp2_ensembl[, "ensembl_gene_sp2"],
      attributes = c("ensembl_gene_id", idtype_sp2),
      mart = dataset_sp2
    ) %>%
      setNames(c("ensembl_gene_id", "id_sp2"))
  }

  # Merge species 1 gene_ids with species 2 gene_ids of type idtype_sp2
  # - removing the intermediate ensembl-ids in species 2
  sp1_to_sp2 <- merge(
    x = sp1_to_sp2_ensembl,
    y = sp2_ensembl_to_output,
    by.x = "ensembl_gene_sp2",
    by.y = "ensembl_gene_id",
    all.x = TRUE
  ) %>%
    dplyr::select_(.dots = c("id_sp1", "id_sp2")) %>%
    setNames(c("id_sp1", "id_sp2")) %>%
    unique()

  # Remove any rows where the sp2.entrez is NA AND the sp1.entrez has at least
  # one non-NA value in sp2.entrez
  has_valid <- with(
    sp1_to_sp2,
    id_sp1[which(!is.na(id_sp2))]
  )

  sp1_to_sp2 %>%
    dplyr::filter_(~id_sp1 %in% has_valid &
      !is.na(id_sp2)) %>%
    dplyr::bind_rows(
      dplyr::filter_(sp1_to_sp2, ~!(id_sp1 %in% has_valid))
    ) %>%
    dplyr::arrange_("id_sp1", "id_sp2") %>%
    as.data.frame(stringsAsFactors = FALSE)
}

###############################################################################

#' map_to_ensembl_homologues_with_biomart
#'
#' A function to map from a set of `gene_ids` in one species (sp1) to the
#'   ensembl-gene id of their homologues in a separate species (sp2). The
#'   mappings are performed using biomart databases.
#'
#' @inheritParams   map_to_homologues_oneway
#'
#' @return       A data.frame with columns id_sp1 and ensembl_gene_sp2.
#'
#' @importFrom   magrittr      %>%
#' @importFrom   dplyr         arrange_
#'
#' @include      homologue_db.R
#'

map_to_ensembl_homologues_with_biomart <- function(gene_ids = character(0),
                                                   dataset_sp1 = NULL,
                                                   sp2 = "mmusculus",
                                                   idtype_sp1 =
                                                     "ensembl_gene_id") {
  results_empty <- data.frame(
    id_sp1 = character(0),
    ensembl_gene_sp2 = character(0),
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
  #   a) be a Mart; b) have a homologue field for sp2; and c) have a field for
  #   the input idtype_sp1

  # TODO:
  # - move mart-validity checks into select_and_filter

  if (is.null(dataset_sp1)) {
    stop("Mart should be provided as `dataset_sp1`")
  }

  stopifnot(
    is_valid_mart(dataset_sp1, idtype_sp1) &&
      is_valid_mart(dataset_sp1, homologue_field)
  )

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
  #   column names id_sp1 and ensembl_gene_sp2
  data.frame(
    id_sp1 = c(
      id_to_ensembl_sp1_sp2[, idtype_sp1],
      missing_inputs
    ),
    ensembl_gene_sp2 = c(
      id_to_ensembl_sp1_sp2[, homologue_field],
      rep(NA, length(missing_inputs))
    ),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange_(.dots = c("id_sp1", "ensembl_gene_sp2"))
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
