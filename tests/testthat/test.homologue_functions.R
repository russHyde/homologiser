###############################################################################

context("Tests for homology-mapping functions")

###############################################################################

test_that(
  "map_to_homologues_oneway", {
    testthat::skip_if_not_installed("mockery")

    gene_ids <- c("1", "2")
    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = c("A", "B")
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues_oneway(
          gene_ids, mock_mart, mock_mart
        ),
        expected = df(
          id_sp1 = gene_ids,
          id_sp2 = bm_hom$mmusculus_homolog_ensembl_gene
        ),
        info = "ensembl-gene to ensembl-gene oneway homology map"
      ),
      .env = "biomaRt"
    )

    # Ensembl to Entrez one-way homology map
    gene_ids <- c("1", "2")

    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = c("A", "B")
    )

    bm_id_sp2 <- df(
      ensembl_gene_id = bm_hom$mmusculus_homolog_ensembl_gene,
      entrezgene = c("entrezA", "entrezB")
    )

    with_mock(
      getBM = mockery::mock(bm_hom, bm_id_sp2),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues_oneway(
          gene_ids, mock_mart, mock_mart,
          idtype_sp2 = "entrezgene"
        ),
        expected = df(
          id_sp1 = gene_ids,
          id_sp2 = c("entrezA", "entrezB")
        ),
        info = "ensembl to entrez-gene oneway homology map"
      ),
      .env = "biomaRt"
    )
  }
)

###############################################################################

test_that(
  "map_to_ensembl_homologues_with_biomart: invalid inputs", {

    # An error should arise from passing in a NULL or non-vector gene list
    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = NULL,
        dataset_sp1 = mock_mart
      ),
      info = paste(
        "NULL gene_ids throw an error in",
        "map_to_ensembl_homologues_with_biomart"
      )
    )

    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = data.frame(gene_ids = c("1000", "1234", "11111")),
        dataset_sp1 = mock_mart
      ),
      info = paste(
        "Non-vector gene_ids throw an error in",
        "map_to_ensembl_homologues_with_biomart"
      )
    )

    # Input gene id type should be either ensembl_gene_id or entrezgene
    id_type_error <- "not_an_id_type"
    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = c("ABC", "123"),
        dataset_sp1 = mock_mart,
        idtype_sp1 = id_type_error
      ),
      info = "invalid id_type in the input"
    )

    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = c("123", "234"),
        dataset_sp1 = "not-a-mart-object"
      ),
      info = "A valid Mart should be passed if dataset_sp1 is non-null"
    )

    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = c("123", "234"),
        dataset_sp1 = NULL
      ),
      info = "dataset_sp1 can not be NULL"
    )

    testthat::skip_if_not_installed("mockery")

    # A valid biomaRt object should be provided for species 1
    # And a homologue column for species 2 should be present in the
    #   biomaRt dataset for species 1

    with_mock(
      listAttributes = mockery::mock(mock_attribs),
      expect_error(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = c("ABC", "123"),
          idtype_sp1 = "entrezgene",
          sp2 = "not_present_in_the_attributes",
          dataset_sp1 = mock_mart
        ),
        info = "Non-species sp2 in map_to_ensembl_homologues_with_biomart"
      ),
      .env = "biomaRt"
    )
  }
)

###############################################################################

test_that(
  "map_to_ensembl_homologues_with_biomart: valid input", {

    # The default return value is a dataframe with two cols and no
    # entries
    expect_empty <- df(
      id_sp1 = character(0),
      ensembl_gene_sp2 = character(0)
    )

    # An empty gene list should result from passing in an empty gene
    # list
    expect_equal(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = character(0)
      ),
      expected = expect_empty,
      info = "Empty gene list should return an empty map"
    )

    expect_equal(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = character(0),
        idtype_sp1 = "entrezgene"
      ),
      expected = expect_empty,
      info = "Empty entrezgene list should return an empty map"
    )

    expect_equal(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = character(0),
        idtype_sp1 = "ensembl_gene_id"
      ),
      expected = expect_empty,
      info = "Empty ensemblgene list should return an empty map"
    )

    testthat::skip_if_not_installed("mockery")

    # Map ensembl gene ids from human to ensembl gene ids in mouse
    # - 1:1 mappings only
    # - results should be alpha-numerically sorted by ID in human
    gene_ids <- c("ENSG00000134294", "ENSG00000117020")

    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000022462", "ENSMUSG00000019699"
      )
    )

    expect <- df(
      id_sp1 = c("ENSG00000117020", "ENSG00000134294"),
      ensembl_gene_sp2 = c("ENSMUSG00000019699", "ENSMUSG00000022462")
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          idtype_sp1 = "ensembl_gene_id",
          sp2 = "mmusculus"
        ),
        expected = expect,
        info = "1:1 mapping of gene homologues - resorting required (hs -> mm)"
      ),
      .env = "biomaRt"
    )

    # map ensembl gene ids from mouse to ensembl gene ids in human
    # - 1:1 mappings only
    # - results should be alpha-numerically sorted by ID in mouse
    gene_ids <- c("ENSMUSG00000019699", "ENSMUSG00000022462")
    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      hsapiens_homolog_ensembl_gene = c("ENSG00000117020", "ENSG00000134294")
    )
    expect <- df(
      id_sp1 = bm_hom$ensembl_gene_id,
      ensembl_gene_sp2 = bm_hom$hsapiens_homolog_ensembl_gene
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          idtype_sp1 = "ensembl_gene_id",
          sp2 = "hsapiens"
        ),
        expected = expect,
        info = paste(
          "1:1 mapping of gene homologues - no resorting required (mm -> hs)"
        )
      ),
      .env = "biomaRt"
    )

    # 1:1 mappings from Entrez IDs in sp1 to Ensembl ID in sp2
    gene_ids <- c("7316", "10000", "54407")

    bm_id <- df(
      ensembl_gene_id = c(
        "ENSG00000150991", "ENSG00000117020", "ENSG00000134294"
      ),
      entrezgene = c("7316", "10000", "54407")
    )
    bm_hom <- df(
      ensembl_gene_id = c(
        "ENSG00000150991", "ENSG00000117020", "ENSG00000134294"
      ),
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000008348", "ENSMUSG00000019699", "ENSMUSG00000022462"
      )
    )
    expect <- df(
      id_sp1 = c("10000", "54407", "7316"),
      ensembl_gene_sp2 = c(
        "ENSMUSG00000019699", "ENSMUSG00000022462", "ENSMUSG00000008348"
      )
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          idtype_sp1 = "entrezgene",
          sp2 = "mmusculus"
        ),
        expected = expect,
        info = "Human entrez ids to mouse ensembl ids (1:1 mapping)"
      ),
      .env = "biomaRt"
    )

    # 1:Many mappings from human to mouse
    gene_ids <- "8813"

    bm_id <- df(
      ensembl_gene_id = "ENSG00000000419",
      entrezgene = "8813"
    )

    bm_hom <- df(
      ensembl_gene_id = rep("ENSG00000000419", 2),
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000078919", "ENSMUSG00000093752"
      )
    )
    expect <- df(
      id_sp1 = c("8813", "8813"),
      ensembl_gene_sp2 = c("ENSMUSG00000078919", "ENSMUSG00000093752")
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          idtype_sp1 = "entrezgene",
          sp2 = "mmusculus"
        ),
        expected = expect,
        info = paste(
          "Human entrez id (8813) that maps to multiple mouse",
          "ensembl-genes"
        )
      ),
      .env = "biomaRt"
    )

    # Many:1 mappings from human to mouse
    # (this is no longer a many:1 in ensembl: 2018-11)
    gene_ids <- c("7543", "7544") # ZFX and ZFY

    bm_id <- df(
      ensembl_gene_id = c("ENSG00000005889", "ENSG00000067646"),
      entrezgene = gene_ids
    )

    bm_hom <- df(
      ensembl_gene_id = bm_id$ensembl_gene_id,
      mmusculus_homolog_ensembl_gene = rep("ENSMUSG00000079509", 2)
    )

    expect <- df(
      id_sp1 = gene_ids,
      ensembl_gene_sp2 = bm_hom$mmusculus_homolog_ensembl_gene
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          idtype_sp1 = "entrezgene",
          sp2 = "mmusculus"
        ),
        expected = expect,
        info = paste(
          "Human entrez ids that map many:1 into mouse ensembl-genes"
        )
      ),
      .env = "biomaRt"
    )

    # Gene that does map from human to mouse
    # - the gene should still be present in the output
    gene_ids <- c("10000", "3", "54407")

    bm_id <- df(
      ensembl_gene_id = c(
        "ENSG00000117020", "ENSG00000256069", "ENSG00000134294"
      ),
      entrezgene = gene_ids
    )

    bm_hom <- df(
      ensembl_gene_id = bm_id$ensembl_gene_id,
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000019699", "", "ENSMUSG00000022462"
      )
    )

    expect <- df(
      id_sp1 = c("10000", "3", "54407"),
      ensembl_gene_sp2 = c(
        "ENSMUSG00000019699", NA, "ENSMUSG00000022462"
      )
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          idtype_sp1 = "entrezgene",
          sp2 = "mmusculus"
        ),
        expected = expect,
        info = paste(
          "Human gene that doesn't map to a mouse gene"
        )
      ),
      .env = "biomaRt"
    )
  }
)

###############################################################################

test_that(
  "map_to_homologues", {
    # Map ensembl to ensembl: human to mouse, single one-to-one example

    testthat::skip_if_not_installed("mockery")

    # No genes
    expect <- df(
      id_sp1 = character(0),
      id_sp2 = character(0)
    )

    with_mock(
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = character(0),
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "hsapiens",
          sp2 = "mmusculus",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id"
        ),
        expected = expect,
        info = "No genes: should return a zero-row dataframe"
      ),
      .env = "biomaRt"
    )

    # Ensembl to ensembl
    gene_ids <- "ENSG00000134294"
    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = "ENSMUSG00000022462"
    )
    expect <- df(
      id_sp1 = "ENSG00000134294",
      id_sp2 = "ENSMUSG00000022462"
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "hsapiens",
          sp2 = "mmusculus",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id"
        ),
        expected = expect,
        info = paste(
          "Ensembl-to-ensembl homologues; human to mouse singly mapping",
          "example"
        )
      ),
      .env = "biomaRt"
    )

    # Then reversed:
    gene_ids <- "ENSMUSG00000022462"
    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      hsapiens_homolog_ensembl_gene = "ENSG00000134294"
    )
    expect <- df(
      id_sp1 = "ENSMUSG00000022462",
      id_sp2 = "ENSG00000134294"
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "mmusculus",
          sp2 = "hsapiens",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id"
        ),
        expected = expect,
        info = paste(
          "Ensembl-to-ensembl homologues; mouse to human singly mapping",
          "example"
        )
      ),
      .env = "biomaRt"
    )

    # Human ensembl gene with no mouse orthologue
    # - note that default mapping is ensembl:ensembl and human:mouse
    # - ENSG00000284192 was found by searching in biomart with filters =
    # "with_mmusculus_homolog", values = FALSE and then keeping only
    # those genes with an entrez gene id
    gene_ids <- "ENSG00000284192"

    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = NA_character_
    )

    expect <- df(
      id_sp1 = gene_ids,
      id_sp2 = NA_character_
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "hsapiens",
          sp2 = "mmusculus",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id"
        ),
        expected = expect,
        info = "Non-mapping human gene, ensembl to ensembl"
      ),
      .env = "biomaRt"
    )

    # Human ensembl gene with multiple mouse orthologues
    gene_ids <- "ENSG00000002726"

    bm_hom <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000029811", "ENSMUSG00000029813", "ENSMUSG00000039215",
        "ENSMUSG00000068536"
      )
    )

    expect <- df(
      id_sp1 = bm_hom$ensembl_gene_id,
      id_sp2 = bm_hom$mmusculus_homolog_ensembl_gene
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "hsapiens",
          sp2 = "mmusculus",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id"
        ),
        expected = expect,
        info = "One-to-many mapping human gene, ensembl-to-ensembl"
      ),
      .env = "biomaRt"
    )
  }
)

###############################################################################

test_that(
  "map_to_homologues: one_to_one", {
    testthat::skip_if_not_installed("mockery")

    # One-to-one ensembl to ensembl: human to mouse, single one-to-one
    # example
    gene_ids <- "ENSG00000134294"
    bm_hom_sp1 <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = "ENSMUSG00000022462"
    )
    bm_hom_sp2 <- df(
      ensembl_gene_id = bm_hom_sp1$mmusculus_homolog_ensembl_gene,
      hsapiens_homolog_ensembl_gene = bm_hom_sp1$ensembl_gene_id
    )
    expect <- df(
      id_sp1 = bm_hom_sp1$ensembl_gene_id,
      id_sp2 = bm_hom_sp1$mmusculus_homolog_ensembl_gene
    )

    with_mock(
      getBM = mockery::mock(bm_hom_sp1, bm_hom_sp2),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "hsapiens",
          sp2 = "mmusculus",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id",
          one_to_one = TRUE
        ),
        expected = expect,
        info = paste(
          "one_to_one ensembl-to-ensembl homologues; singly mapping example"
        )
      ),
      .env = "biomaRt"
    )

    # 1:many
    gene_ids <- "ENSG00000002726"

    bm_hom_sp1 <- df(
      ensembl_gene_id = gene_ids,
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000029811", "ENSMUSG00000029813", "ENSMUSG00000039215",
        "ENSMUSG00000068536"
      )
    )

    bm_hom_sp2 <- df(
      ensembl_gene_id = bm_hom_sp1$mmusculus_homolog_ensembl_gene,
      hsapiens_homolog_ensembl_gene = bm_hom_sp1$ensembl_gene_id
    )

    expect <- df(
      id_sp1 = bm_hom_sp1$ensembl_gene_id[1],
      id_sp2 = NA_character_
    )

    with_mock(
      getBM = mockery::mock(bm_hom_sp1, bm_hom_sp2),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "hsapiens",
          sp2 = "mmusculus",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id",
          one_to_one = TRUE
        ),
        expected = expect,
        info = paste(
          "one_to_one ensembl-to-ensembl homologues; 1:many example"
        )
      ),
      .env = "biomaRt"
    )

    # Many to one example; note *29811 is one of four mouse homologues
    # of ENSG0000002726

    gene_ids <- "ENSMUSG00000029811"

    bm_hom_sp1 <- df(
      ensembl_gene_id = gene_ids,
      hsapiens_homolog_ensembl_gene = "ENSG00000002726"
    )

    bm_hom_sp2 <- df(
      ensembl_gene_id = bm_hom_sp1$hsapiens_homolog_ensembl_gene,
      mmusculus_homolog_ensembl_gene = c(
        "ENSMUSG00000029811", "ENSMUSG00000029813", "ENSMUSG00000039215",
        "ENSMUSG00000068536"
      )
    )

    expect <- df(
      id_sp1 = bm_hom_sp1$ensembl_gene_id[1],
      id_sp2 = NA_character_
    )

    with_mock(
      getBM = mockery::mock(bm_hom_sp1, bm_hom_sp2),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_homologues(
          gene_ids = gene_ids,
          dataset_sp1 = mock_mart,
          dataset_sp2 = mock_mart,
          sp1 = "mmusculus",
          sp2 = "hsapiens",
          idtype_sp1 = "ensembl_gene_id",
          idtype_sp2 = "ensembl_gene_id",
          one_to_one = TRUE
        ),
        expected = expect,
        info = paste(
          "one_to_one ensembl-to-ensembl homologues; 1:many example"
        )
      ),
      .env = "biomaRt"
    )
  }
)

###############################################################################

test_that("keep_complete_biomart_results", {
  # TODO: refactor select_and_filter - put the empty -> NA conversion and
  # filtering steps in keep_complete_biomart_results
})

###############################################################################

test_that("select_and_filter", {
  testthat::skip_if_not_installed("mockery")

  my_bm <- df(
    ensembl_gene_id = "ENSG01234567890",
    entrezgene = "12345"
  )

  my_missingness_bm <- df(
    ensembl_gene_id = c("ENSG01234567890", "ENSG09876543210"),
    entrezgene = c("12345", "")
  )

  with_mock(
    getBM = mockery::mock(my_bm),
    expect_equal(
      select_and_filter(
        "ensembl_gene_id",
        values = NA_character_, colnames(my_bm), "mock_mart"
      ),
      df(
        ensembl_gene_id = character(0), entrezgene = character(0)
      ),
      info = "no valid values in the input: output is empty"
    ),
    .env = "biomaRt"
  )

  with_mock(
    getBM = mockery::mock(my_bm),
    expect_equal(
      select_and_filter(
        "ensembl_gene_id", my_bm$ensembl_gene_id, colnames(my_bm), "mock_mart"
      ),
      my_bm,
      info = "all returned values are non-empty: no filtering is needed"
    ),
    .env = "biomaRt"
  )

  with_mock(
    getBM = mockery::mock(my_missingness_bm),
    expect_equivalent(
      select_and_filter(
        "ensembl_gene_id", my_missingness_bm$ensembl_gene_id,
        colnames(my_missingness_bm), "mock_mart"
      ),
      my_missingness_bm[1, 1:2, drop = FALSE],
      info = paste(
        "`select_and_filter` should remove any incomplete rows from the",
        "biomart data"
      )
    ),
    .env = "biomaRt"
  )
})

###############################################################################
