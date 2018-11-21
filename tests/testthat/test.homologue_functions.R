###############################################################################

context("Tests for homology-mapping functions")

###############################################################################
## helpers
df <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}

###############################################################################

# TODO: mock out all relevant tests

test_that(
  "use_default_mart", {

    #      # Test that a valid Mart object is returned
    #      expect_true(
    #        is(use_default_mart(), "Mart"),
    #        info = "test that a valid mart is returned"
    #      )

    #      # Test that a valid Mart object is returned - using named args
    #      expect_true(
    #        is(
    #          use_default_mart(
    #            sp = "hsapiens",
    #            host = "www.ensembl.org",
    #            mart_name = "ENSEMBL_MART_ENSEMBL"
    #          ),
    #          "Mart"
    #        ),
    #        info = "test that a valid mart is returned"
    #      )

    #      # species name, hostname and martname should all be valid in
    #      # use_default_mart
    #      expect_error(
    #        use_default_mart(
    #          sp = "NOT_A_SPECIES",
    #          host = "www.ensembl.org",
    #          mart_name = "ENSEMBL_MART_ENSEMBL"
    #        ),
    #        info = "non-standard species name"
    #      )
    #      expect_error(
    #        use_default_mart(
    #          sp = "hsapiens",
    #          host = "NOT_A_URL",
    #          mart_name = "ENSEMBL_MART_ENSEMBL"
    #        ),
    #        info = "non-standard host URL"
    #      )
    #      expect_error(
    #        use_default_mart(
    #          sp = "hsapiens",
    #          host = "www.ensembl.org",
    #          mart_name = "NOT_A_MARTNAME"
    #        ),
    #        info = "non-standard mart_name"
    #      )
  }
)
###############################################################################


test_that(
  "map_to_ensembl_homologues_with_biomart: invalid inputs", {

    # An error should arise from passing in a NULL or non-vector gene list
    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = NULL
      ),
      info = paste(
        "NULL gene_ids throw an error in",
        "map_to_ensembl_homologues_with_biomart"
      )
    )

    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = data.frame(gene_ids = c("1000", "1234", "11111"))
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
        dataset_sp1 = NULL,
        sp1 = NULL
      ),
      info = "A way to define a valid Mart for species 1 should be provided"
    )

    # Mock biomaRt objects for use in unit-testing
    mock_mart <- structure(.Data = list(), class = "Mart")

    mock_attribs <- df(
      name = c(
        "ensembl_gene_id", "entrezgene", "mmusculus_homolog_ensembl_gene",
        "hsapiens_homolog_ensembl_gene"
      ),
      descriptions = c(
        "Gene ID", "External Gene ID", "Mouse homologues", "Human homologues"
      )
    )
    # A valid biomaRt object should be constructible for species 1
    # And a homologue column for species 2 should be present in the
    #   biomaRt dataset for species 1

    with_mock(
      listAttributes = mockery::mock(mock_attribs),
      expect_error(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = c("ABC", "123"),
          idtype_sp1 = "entrezgene",
          sp1 = "mmusculus",
          sp2 = "not_present_in_the_attributes"
        ),
        info = "Non-species sp2 in map_to_ensembl_homologues_with_biomart"
      ),
      .env = "biomaRt"
    )

    # TODO: all biomart calls should be mocked out but I can't work out how
    # to mock out the use_default_mart calls
    # - There should be no calls to use_default_mart within
    # map_to_ensembl_homologues_with_biomart; this is a non-exported function
    # and so we can rewrite it to always receive a Mart

    # # input mart/host/species1 should imply a valid biomart dataset
    # expect_error(
    #   object = map_to_ensembl_homologues_with_biomart(
    #     gene_ids = c("ABC", "123"),
    #     idtype_sp1 = "entrezgene",
    #     sp1 = "I'm not a species"
    #   ),
    #   info = "Non-species sp1 in map_to_ensembl_homologues_with_biomart"
    # )

    with_mock(
      url.exists = mockery::mock(FALSE),
      expect_error(
        object = map_to_ensembl_homologues_with_biomart(
         gene_ids = c("ABC", "123"),
         idtype_sp1 = "ensembl_gene_id",
         host = "NOT_A_URL"
       ),
       info = "host-URL is invalid"
      ),
      .env = "RCurl"
    )

    # expect_error(
    #   object = map_to_ensembl_homologues_with_biomart(
    #     gene_ids = c("ABC", "123"),
    #     idtype_sp1 = "ensembl_gene_id",
    #     mart_name = "NOT_A_MARTNAME"
    #   ),
    #   info = "mart_name is invalid"
    # )
  }
)
###############################################################################

test_that(
  "map_to_ensembl_homologues_with_biomart: valid input", {

    # The default return value is a dataframe with two cols and no
    # entries
    expect_empty <- df(
      ID.sp1 = character(0),
      ENSEMBLGENE.sp2 = character(0)
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

    # Mock biomaRt objects for use in unit-testing
    mock_mart <- structure(.Data = list(), class = "Mart")

    mock_attribs <- df(
      name = c(
        "ensembl_gene_id", "entrezgene", "mmusculus_homolog_ensembl_gene",
        "hsapiens_homolog_ensembl_gene"
      ),
      descriptions = c(
        "Gene ID", "External Gene ID", "Mouse homologues", "Human homologues"
      )
    )

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
      ID.sp1 = c("ENSG00000117020", "ENSG00000134294"),
      ENSEMBLGENE.sp2 = c("ENSMUSG00000019699", "ENSMUSG00000022462")
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          idtype_sp1 = "ensembl_gene_id",
          sp1 = "hsapiens",
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
      ID.sp1 = bm_hom$ensembl_gene_id,
      ENSEMBLGENE.sp2 = bm_hom$hsapiens_homolog_ensembl_gene
    )

    with_mock(
      getBM = mockery::mock(bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          idtype_sp1 = "ensembl_gene_id",
          sp1 = "mmusculus",
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
      ID.sp1 = c("10000", "54407", "7316"),
      ENSEMBLGENE.sp2 = c(
        "ENSMUSG00000019699", "ENSMUSG00000022462", "ENSMUSG00000008348"
      )
    )

    expect_equal(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = gene_ids,
        idtype_sp1 = "entrezgene",
        sp1 = "hsapiens",
        sp2 = "mmusculus"
      ),
      expected = expect,
      info = "Human entrez ids to mouse ensembl ids (1:1 mapping)"
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          idtype_sp1 = "entrezgene",
          sp1 = "hsapiens",
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
      ID.sp1 = c("8813", "8813"),
      ENSEMBLGENE.sp2 = c("ENSMUSG00000078919", "ENSMUSG00000093752")
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          idtype_sp1 = "entrezgene",
          sp1 = "hsapiens",
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
      ID.sp1 = gene_ids,
      ENSEMBLGENE.sp2 = bm_hom$mmusculus_homolog_ensembl_gene
    )

    with_mock(
      getBM = mockery::mock(bm_id, bm_hom),
      useMart = mockery::mock(mock_mart, cycle = TRUE),
      listAttributes = mockery::mock(mock_attribs, cycle = TRUE),
      expect_equal(
        object = map_to_ensembl_homologues_with_biomart(
          gene_ids = gene_ids,
          idtype_sp1 = "entrezgene",
          sp1 = "hsapiens",
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
      ID.sp1 = c("10000", "3", "54407"),
      ENSEMBLGENE.sp2 = c(
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
          idtype_sp1 = "entrezgene",
          sp1 = "hsapiens",
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
  "map_to_homologues_with_biomart", {

    # TODO: rename all test-variables in snake_case
    # expect the results to look like
    # data.frame(ENTREZID.sp1 = character.vec, ENTREZID.sp2 = character.vec)
    # Empty gene list should return an empty map:
    gene_ids_empty <- character(0)
    # expect.empty <- data.frame(
    #  ENTREZID.sp1 = character(0),
    #  ENTREZID.sp2 = character(0),
    #  stringsAsFactors = FALSE
    # )
    #      result.empty <- map_to_homologues_with_biomart(gene_ids_empty)
    #      expect_equal(
    #        object = expect.empty,
    #        expected = result.empty,
    #        info = "Empty gene list should return an empty map"
    #      )

    #      # Null gene list should return an empty map
    #      gene_ids_null <- NULL
    #      expect.null <- data.frame(
    #        ENTREZID.sp1 = character(0),
    #        ENTREZID.sp2 = character(0),
    #        stringsAsFactors = FALSE
    #      )
    #      result.null <- map_to_homologues_with_biomart(gene_ids_null)
    #      expect_equal(
    #        object = expect.null,
    #        expected = result.null,
    #        info = "Null gene list should return an empty map"
    #      )

    #      # Gene.ids should be character,
    #      #   so numeric or factor input should die
    #      gene_ids_numeric <- 1000
    #      gene_ids_factor <- as.factor("1000")
    #      expect_error(
    #        map_to_homologues_with_biomart(gene_ids_numeric),
    #        info = "Gene ids should be a character vector: Numeric input"
    #      )
    #      expect_error(
    #        map_to_homologues_with_biomart(gene_ids_factor),
    #        info = "Gene ids should be a character vector: Factor input"
    #      )

    #      # A human gene that doesn"t have a corresponding ensembl-mouse id
    #      #   cannot be mapped between species using ensembl/biomart
    #      #   and should return NA
    #      # Entrez gene "3" maps to "NA" ensembl prot:
    #      gene_ids_without.mouse <- c("3")
    #      expect_without_mouse <- data.frame(
    #        ENTREZID.sp1 = gene_ids_without.mouse,
    #        ENTREZID.sp2 = as.character(NA),
    #        stringsAsFactors = FALSE
    #      )
    #      expect_equal(
    #        object = map_to_homologues_with_biomart(gene_ids_without.mouse),
    #        expected = expect_without_mouse,
    #        info = "Human gene without ensembl-mouse id"
    #      )

    #      # A human gene that has a mouse orthologue
    #      # Entrez:1
    #      # maps to ENSG00000121410
    #      # maps to mouse ENSMUSG00000022347 and ENSMUSG00000112930
    #      # map to        Entrez:117586      and Entrez:NA
    #      # If a gene maps to some genuine orthologue IDs and to some NAs, the
    #      # NAs should not be present in the output
    #      expect_equal(
    #        object = map_to_homologues_with_biomart("1"),
    #        expected = data.frame(
    #          ENTREZID.sp1 = "1",
    #          ENTREZID.sp2 = "117586",
    #          stringsAsFactors = FALSE
    #        ),
    #        info = paste(
    #          "Human gene that maps to both valid-ID and NA in mouse, the NA",
    #          "should be dropped when a valid ID is available"
    #          )
    #      )

    #      # Human genes that do have mouse orthologues:
    #      # nb, "2" and "41" did not map to an orthologue using inparanoid
    #      #     though Mm orthologues exist on ensembl.org
    #      gene_ids_with.mouse <- c("1", "151", "2", "33", "41")
    #      expect_with_mouse <- data.frame(
    #        ENTREZID.sp1 = gene_ids_with.mouse,
    #        ENTREZID.sp2 = c("117586", "11552", "232345", "11363", "11419"),
    #        stringsAsFactors = FALSE
    #      )
    #      expect_equal(
    #        object = map_to_homologues_with_biomart(gene_ids_with.mouse),
    #        expected = expect_with_mouse,
    #        info = "Human genes with mouse orthologues"
    #      )

    #      # Human gene that maps to multiple mouse orthologues
    #      # note, additional orthologue, compared to inparanoid searches
    #      # note, 102641229 is not suggested as an orthologue on
    #      #   www.ensembl.org but does return under biomaRt queries (?
    #      #   predicted assembly
    #      # THIS TEST FAILS AS OF 14/12/2015
    #      # - Wrote up a new test of this type of gene using "26" as input
    #      # - See below
    #      # gene_ids_with_multi_mouse <- "8363"
    #      # expect_with_multi.mouse <- data.frame(
    #      #  ENTREZID.sp1 = "8363",
    #      #  ENTREZID.sp2 = c("102641229", "319155", "326620"),
    #      #  stringsAsFactors = FALSE
    #      #  )

    #      # Human gene that maps to multiple mouse ensembl gene ids via a
    #      # single human ensembl gene id
    #      gene_ids_with_multi_mouse <- "26"
    #      expect_with_multi.mouse <- data.frame(
    #        ENTREZID.sp1 = rep("26", 4),
    #        ENTREZID.sp2 = c("243376", "243377", "69761", "76507"),
    #        stringsAsFactors = FALSE
    #      )
    #      expect_equal(
    #        object = map_to_homologues_with_biomart(
    #          gene_ids_with_multi_mouse
    #        ),
    #        expected = expect_with_multi.mouse,
    #        info = paste(
    #          "Human entrez id that maps to multiple mouse entrez ids",
    #          "(via a single HS-ensembl and multiple MM-ensembl ids)"
    #        )
    #      )

    #      # Human genes that map to the same mouse orthologues
    #      # "7314"/UBB and "7316"/UBC map to both "22187"/Ubb and "22190"/Ubc
    #      #
    #      # The test that uses UBB/UBC no longer works as of 14/12/2015
    #      #   since 7314 maps to a single gene and 7316 maps to a single gene
    #      #   now gene_ids_many_to_many <- c("7314", "7316")
    #      # expect.many_to_many <- data.frame(
    #      #   ENTREZID.sp1 = rep(gene_ids_many_to_many, each = 2),
    #      #   ENTREZID.sp2 = rep(c("22187", "22190"), 2),
    #      #   stringsAsFactors = FALSE
    #      #   )
    #      # Hence, I rewrote the method using the human alkaline phosphatases
    #      # ALPP and ALPPL2
    #      # Each of these genes maps to 3 different mouse gene ids (in biomart
    #      # on 14/12/15)
    #      gene_ids_many_to_many <- c("250", "251")
    #      expect.many_to_many <- data.frame(
    #        ENTREZID.sp1 = rep(c("250", "251"), each = 3),
    #        ENTREZID.sp2 = rep(c("11648", "11650", "76768"), times = 2),
    #        stringsAsFactors = FALSE
    #      )
    #      expect_equal(
    #        object = map_to_homologues_with_biomart(
    #          gene_ids_many_to_many
    #        ),
    #        expected = expect.many_to_many,
    #        info = "Genes that map to the same multiple genes"
    #      )
  }
)

###############################################################################

test_that(
  "map_to_homologues", {
    #
    #      # Map ensembl to ensembl: human to mouse, single one-to-one example
    #      expect_equal(
    #        object = map_to_homologues(
    #          gene_ids = "ENSG00000134294",
    #          sp1 = "hsapiens",
    #          sp2 = "mmusculus",
    #          idtype_sp1 = "ensembl_gene_id",
    #          idtype_sp2 = "ensembl_gene_id"
    #        ),
    #        expected = data.frame(
    #          ID.sp1 = "ENSG00000134294",
    #          ID.sp2 = "ENSMUSG00000022462",
    #          stringsAsFactors = FALSE
    #        ),
    #        info = paste(
    #          "Ensembl-to-ensembl homologues; human to mouse singly",
    #          "mapping example"
    #        )
    #      )
    #      # Then reversed:
    #      expect_equal(
    #        object = map_to_homologues(
    #          gene_ids = "ENSMUSG00000022462",
    #          sp1 = "mmusculus",
    #          sp2 = "hsapiens",
    #          idtype_sp1 = "ensembl_gene_id",
    #          idtype_sp2 = "ensembl_gene_id"
    #        ),
    #        expected = data.frame(
    #          ID.sp1 = "ENSMUSG00000022462",
    #          ID.sp2 = "ENSG00000134294",
    #          stringsAsFactors = FALSE
    #        ),
    #        info = paste(
    #          "Ensembl-to-ensembl homologues; mouse to human singly",
    #          "mapping example"
    #        )
    #      )
    #      # Human ensembl gene with no mouse orthologue
    #      # - note that default mapping is ensembl:ensembl and human:mouse
    #      # - ENSG00000284192 was found by searching in biomart with filters =
    #      # "with_mmusculus_homolog", values = FALSE and then keeping only
    #      # those genes with an entrez gene id
    #      expect_equal(
    #        object = map_to_homologues("ENSG00000284192"),
    #        expected = data.frame(
    #          ID.sp1 = "ENSG00000284192",
    #          ID.sp2 = as.character(NA),
    #          stringsAsFactors = FALSE
    #        ),
    #        info = "Non-mapping human gene, ensembl to ensembl"
    #      )
    #      # Human ensembl gene with multiple mouse orthologues
    #      expect_equal(
    #        object = map_to_homologues("ENSG00000002726"),
    #        expected = data.frame(
    #          ID.sp1 = "ENSG00000002726",
    #          ID.sp2 = c(
    #            "ENSMUSG00000029811",
    #            "ENSMUSG00000029813",
    #            "ENSMUSG00000039215",
    #            "ENSMUSG00000068536"
    #          ),
    #          stringsAsFactors = FALSE
    #        ),
    #        info = "One-to-many mapping human gene, ensembl-to-ensembl"
    #      )

    #      # Would like a gene that maps many:many for use in testing
    #      # one.to.one = TRUE
  }
)

###############################################################################

test_that(
  "map_to_homologues: one_to_one", {
    #      # One-to-one ensembl to ensembl: human to mouse, single one-to-one
    #      # example
    #      expect_equal(
    #        object = map_to_homologues(
    #          gene_ids = "ENSG00000134294",
    #          one_to_one = TRUE
    #        ),
    #        expected = data.frame(
    #          ID.sp1 = "ENSG00000134294",
    #          ID.sp2 = "ENSMUSG00000022462",
    #          stringsAsFactors = FALSE
    #        ),
    #        info = paste(
    #          "one_to_one ensembl-to-ensembl homologues; singly mapping",
    #          "example"
    #        )
    #      )
    #      expect_equal(
    #        object = map_to_homologues(
    #          gene_ids = "ENSG00000002726",
    #          one_to_one = TRUE
    #        ),
    #        expected = DF(
    #          ID.sp1 = "ENSG00000002726",
    #          ID.sp2 = as.character(NA)
    #        ),
    #        info = "one_to_one ensembl-to-ensembl homologues; 1:many example"
    #      )
    #      expect_equal(
    #        # Many to one example; note *29811 is one of four mouse homologues
    #        # of ENSG0000002726
    #        object = map_to_homologues(
    #          gene_ids = "ENSMUSG00000029811",
    #          sp1 = "mmusculus",
    #          sp2 = "hsapiens",
    #          one_to_one = TRUE
    #        ),
    #        expected = DF(
    #          ID.sp1 = "ENSMUSG00000029811",
    #          ID.sp2 = as.character(NA)
    #        ),
    #        info = "one_to_one ensembl-to-ensembl homologues; many:1 example"
    #      )
  }
)

###############################################################################

test_that("get_ensembl_homologue_field", {
  expect_equal(
    get_ensembl_homologue_field("hsapiens"),
    "hsapiens_homolog_ensembl_gene",
    info = "correct formatting of homology field for a given species"
  )

  expect_error(
    get_ensembl_homologue_field(),
    info = "homologue-field input should be defined"
  )

  expect_error(
    get_ensembl_homologue_field(NULL),
    info = "homologue-field input should be non-null"
  )

  expect_error(
    get_ensembl_homologue_field(123),
    info = "homologue-field is only defined for strings"
  )

  expect_error(
    get_ensembl_homologue_field("Homo sapiens"),
    info = "no spaces allowed in the homologue-field"
  )
})

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

test_that("keep_complete_biomart_results", {
  # TODO: refactor select_and_filter - put the empty -> NA conversion and
  # filtering steps in keep_complete_biomart_results
})

test_that("select_and_filter", {
  testthat::skip_if_not(requireNamespace("mockery", quietly = TRUE))

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

# TODO: mock out all relevant tests

test_that(
  "use_default_mart", {

    #      # Test that a valid Mart object is returned
    #      expect_true(
    #        is(use_default_mart(), "Mart"),
    #        info = "test that a valid mart is returned"
    #      )

    #      # Test that a valid Mart object is returned - using named args
    #      expect_true(
    #        is(
    #          use_default_mart(
    #            sp = "hsapiens",
    #            host = "www.ensembl.org",
    #            mart_name = "ENSEMBL_MART_ENSEMBL"
    #          ),
    #          "Mart"
    #        ),
    #        info = "test that a valid mart is returned"
    #      )

    #      # species name, hostname and martname should all be valid in
    #      # use_default_mart
    #      expect_error(
    #        use_default_mart(
    #          sp = "NOT_A_SPECIES",
    #          host = "www.ensembl.org",
    #          mart_name = "ENSEMBL_MART_ENSEMBL"
    #        ),
    #        info = "non-standard species name"
    #      )
    #      expect_error(
    #        use_default_mart(
    #          sp = "hsapiens",
    #          host = "NOT_A_URL",
    #          mart_name = "ENSEMBL_MART_ENSEMBL"
    #        ),
    #        info = "non-standard host URL"
    #      )
    #      expect_error(
    #        use_default_mart(
    #          sp = "hsapiens",
    #          host = "www.ensembl.org",
    #          mart_name = "NOT_A_MARTNAME"
    #        ),
    #        info = "non-standard mart_name"
    #      )
  }
)
###############################################################################

# TODO: all biomart calls should be mocked out

test_that(
  "map_to_ensembl_homologues_with_biomart: invalid inputs", {

    # An error should arise from passing in a NULL or non-vector gene list
    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = NULL
      ),
      info = paste(
        "NULL gene_ids throw an error in",
        "map_to_ensembl_homologues_with_biomart"
      )
    )

    expect_error(
      object = map_to_ensembl_homologues_with_biomart(
        gene_ids = data.frame(gene_ids = c("1000", "1234", "11111"))
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
        idtype_sp1 = id_type_error
      ),
      info = "invalid id_type in the input"
    )

    #  # A valid biomaRt object should be constructible for species 1
    #  # And a homologue column for species 2 should be present in the biomaRt
    #  #   dataset for species 1
    #  expect_error(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = c("ABC", "123"),
    #      idtype_sp1 = "entrezgene",
    #      sp1 = "mmusculus",
    #      sp2 = "Ri diculous"
    #    ),
    #    info = "Non-species sp2 in map_to_ensembl_homologues_with_biomart"
    #  )

    #  # input mart/host/species1 should imply a valid biomart dataset
    #  expect_error(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = c("ABC", "123"),
    #      idtype_sp1 = "entrezgene",
    #      sp1 = "I'm not a species"
    #    ),
    #    info = "Non-species sp1 in map_to_ensembl_homologues_with_biomart"
    #  )

    #  expect_error(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = c("ABC", "123"),
    #      idtype_sp1 = "ensembl_gene_id",
    #      host = "NOT_A_URL"
    #    ),
    #    info = "host-URL is invalid"
    #  )

    #  expect_error(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = c("ABC", "123"),
    #      idtype_sp1 = "ensembl_gene_id",
    #      mart_name = "NOT_A_MARTNAME"
    #    ),
    #    info = "mart_name is invalid"
    #  )
  }
)

###############################################################################

test_that(
  "map_to_ensembl_homologues_with_biomart: valid input", {

    #  # The default return value is a dataframe with two cols and no entries
    #  expect.empty <- data.frame(
    #    ID.sp1 = character(0),
    #    ENSEMBLGENE.sp2 = character(0),
    #    stringsAsFactors = FALSE
    #  )

    #  # An empty gene list should result from passing in an empty gene list
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = character(0)
    #    ),
    #    expected = expect.empty,
    #    info = "Empty gene list should return an empty map"
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = character(0),
    #      idtype_sp1 = "entrezgene"
    #    ),
    #    expected = expect.empty,
    #    info = "Empty entrezgene list should return an empty map"
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = character(0),
    #      idtype_sp1 = "ensembl_gene_id"
    #    ),
    #    expected = expect.empty,
    #    info = "Empty ensemblgene list should return an empty map"
    #  )

    #  # Map ensembl gene ids from human to ensembl gene ids in mouse
    #  # - 1:1 mappings only
    #  # - results should be alpha-numerically sorted by ID in human
    #  gene_ids <- c(
    #    "ENSG00000134294",
    #    "ENSG00000117020"
    #  )
    #  expect <- data.frame(
    #    ID.sp1 = c("ENSG00000117020", "ENSG00000134294"),
    #    ENSEMBLGENE.sp2 = c("ENSMUSG00000019699", "ENSMUSG00000022462"),
    #    stringsAsFactors = FALSE
    #  )
    #  result <- expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = gene_ids,
    #      idtype_sp1 = "ensembl_gene_id",
    #      sp1 = "hsapiens",
    #      sp2 = "mmusculus"
    #    ),
    #    expected = expect,
    #    info = "1:1 mapping of human ensembl genes to mouse ensembl genes"
    #  )

    #  # map ensembl gene ids from mouse to ensembl gene ids in human
    #  # - 1:1 mappings only
    #  # - results should be alpha-numerically sorted by ID in mouse
    #  gene_ids <- c(
    #    "ENSMUSG00000019699",
    #    "ENSMUSG00000022462"
    #  )
    #  expect <- data.frame(
    #    ID.sp1 = c("ENSMUSG00000019699", "ENSMUSG00000022462"),
    #    ENSEMBLGENE.sp2 = c("ENSG00000117020", "ENSG00000134294"),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = gene_ids,
    #      idtype_sp1 = "ensembl_gene_id",
    #      sp1 = "mmusculus",
    #      sp2 = "hsapiens"
    #    ),
    #    expected = expect,
    #    info = "1:1 mapping of mouse ensembl genes to human ensembl genes"
    #  )

    #  # Entrez IDs in sp1 to Ensembl ID in sp2
    #  # Note that the results should be ordered by ID.sp1 (alphanumerically)
    #  # Note that mappings are volatile; 7316 used to map to
    #  #   ENSMUSG00000019505 but now maps to ENSMUSG00000008348
    #  gene_ids <- c(
    #    "7316", # UBC
    #    "10000", # AKT3
    #    "54407" # SLC38A2
    #  )
    #  expect <- data.frame(
    #    ID.sp1 = c(
    #      "10000",
    #      "54407",
    #      "7316"
    #    ),
    #    ENSEMBLGENE.sp2 = c(
    #      "ENSMUSG00000019699",
    #      "ENSMUSG00000022462",
    #      "ENSMUSG00000008348"
    #    ),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = gene_ids,
    #      idtype_sp1 = "entrezgene",
    #      sp1 = "hsapiens",
    #      sp2 = "mmusculus"
    #    ),
    #    expected = expect,
    #    info = "Human entrez ids to mouse ensembl ids (1:1 mapping ?)"
    #  )

    #  # Check the return from HS:8813 (maps to two mouse-ensembls)
    #  # Note that the ensembl genes should be ordered alphabetically
    #  expect <- data.frame(
    #    ID.sp1 = c("8813", "8813"),
    #    ENSEMBLGENE.sp2 = c("ENSMUSG00000078919", "ENSMUSG00000093752"),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = "8813",
    #      idtype_sp1 = "entrezgene",
    #      sp1 = "hsapiens",
    #      sp2 = "mmusculus"
    #    ),
    #    expected = expect,
    #    info = paste(
    #      "Human entrez id (8813) that maps to multiple mouse ensembl-genes"
    #    )
    #  )

    #  # Check results with HS:7543 and HS:7544; both map to the same mouse
    #  # gene
    #  expect <- data.frame(
    #    ID.sp1 = rep(c("7543", "7544"), each = 3),
    #    ENSEMBLGENE.sp2 = rep(c(
    #      "ENSMUSG00000000103", "ENSMUSG00000053211",
    #      "ENSMUSG00000079509"
    #    ), times = 2),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = c("7543", "7544"),
    #      idtype_sp1 = "entrezgene",
    #      sp1 = "hsapiens",
    #      sp2 = "mmusculus"
    #    ),
    #    expected = expect,
    #    info = paste(
    #      "Human entrez ids (7543/7544) that map to multiple mouse",
    #      "ensembl-genes"
    #    )
    #  )

    #  # Check that a row is present for entrezgene "3" in the output, this
    #  # does not map to an ensembl id in mouse; it should come out sorted
    #  # alphabetically, regardless of not mapping
    #  expect <- data.frame(
    #    ID.sp1 = c("10000", "3", "54407"),
    #    ENSEMBLGENE.sp2 = c(
    #      "ENSMUSG00000019699",
    #      NA,
    #      "ENSMUSG00000022462"
    #    ),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_ensembl_homologues_with_biomart(
    #      gene_ids = c("54407", "3", "10000"),
    #      idtype_sp1 = "entrezgene",
    #      sp1 = "hsapiens",
    #      sp2 = "mmusculus"
    #    ),
    #    expected = expect,
    #    info = paste(
    #      "Human entrez id '3' should not map into mouse (but return NA)",
    #      " - but should still be ordered alphabetically"
    #    )
    #  )
  }
)

###############################################################################

test_that(
  "map_to_homologues_with_biomart", {

    #  # expect the results to look like the results from ..with_inparanoid
    #  # ie, data.frame(ENTREZID.sp1 = character.vec, ENTREZID.sp2 =
    #  # character.vec)
    #  # Empty gene list should return an empty map:
    #  gene_ids_empty <- character(0)
    #  expect.empty <- data.frame(
    #    ENTREZID.sp1 = character(0),
    #    ENTREZID.sp2 = character(0),
    #    stringsAsFactors = FALSE
    #  )
    #  result.empty <- map_to_homologues_with_biomart(gene_ids_empty)
    #  expect_equal(
    #    object = expect.empty,
    #    expected = result.empty,
    #    info = "Empty gene list should return an empty map"
    #  )

    #  # Null gene list should return an empty map
    #  gene_ids_null <- NULL
    #  expect.null <- data.frame(
    #    ENTREZID.sp1 = character(0),
    #    ENTREZID.sp2 = character(0),
    #    stringsAsFactors = FALSE
    #  )
    #  result.null <- map_to_homologues_with_biomart(gene_ids_null)
    #  expect_equal(
    #    object = expect.null,
    #    expected = result.null,
    #    info = "Null gene list should return an empty map"
    #  )

    #  # Gene.ids should be character,
    #  #   so numeric or factor input should die
    #  gene_ids_numeric <- 1000
    #  gene_ids_factor <- as.factor("1000")
    #  expect_error(
    #    map_to_homologues_with_biomart(gene_ids_numeric),
    #    info = "Gene ids should be a character vector: Numeric input"
    #  )
    #  expect_error(
    #    map_to_homologues_with_biomart(gene_ids_factor),
    #    info = "Gene ids should be a character vector: Factor input"
    #  )

    #  # A human gene that doesn"t have a corresponding ensembl-mouse id
    #  #   cannot be mapped between species using ensembl/biomart
    #  #   and should return NA
    #  # Entrez gene "3" maps to "NA" ensembl prot:
    #  gene_ids_without.mouse <- c("3")
    #  expect_without_mouse <- data.frame(
    #    ENTREZID.sp1 = gene_ids_without.mouse,
    #    ENTREZID.sp2 = as.character(NA),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_homologues_with_biomart(gene_ids_without.mouse),
    #    expected = expect_without_mouse,
    #    info = "Human gene without ensembl-mouse id"
    #  )

    #  # A human gene that has a mouse orthologue
    #  # Entrez:1
    #  # maps to ENSG00000121410
    #  # maps to mouse ENSMUSG00000022347 and ENSMUSG00000112930
    #  # map to        Entrez:117586      and Entrez:NA
    #  # If a gene maps to some genuine orthologue IDs and to some NAs, the NAs
    #  # should not be present in the output
    #  expect_equal(
    #    object = map_to_homologues_with_biomart("1"),
    #    expected = data.frame(
    #      ENTREZID.sp1 = "1",
    #      ENTREZID.sp2 = "117586",
    #      stringsAsFactors = FALSE
    #    ),
    #    info = "Human gene that maps to both valid-ID and NA in mouse, the NA
    #  should be dropped when a valid ID is available"
    #  )

    #  # Human genes that do have mouse orthologues:
    #  # nb, "2" and "41" did not map to an orthologue using inparanoid
    #  #     though Mm orthologues exist on ensembl.org
    #  gene_ids_with.mouse <- c("1", "151", "2", "33", "41")
    #  expect_with_mouse <- data.frame(
    #    ENTREZID.sp1 = gene_ids_with.mouse,
    #    ENTREZID.sp2 = c("117586", "11552", "232345", "11363", "11419"),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_homologues_with_biomart(gene_ids_with.mouse),
    #    expected = expect_with_mouse,
    #    info = "Human genes with mouse orthologues"
    #  )

    #  # Human gene that maps to multiple mouse orthologues
    #  # note, additional orthologue, compared to inparanoid searches
    #  # note, 102641229 is not suggested as an orthologue on www.ensembl.org
    #  #     but does return under biomaRt queries (? predicted assembly
    #  # THIS TEST FAILS AS OF 14/12/2015
    #  # - Wrote up a new test of this type of gene using "26" as input
    #  # - See below
    #  # gene_ids_with_multi_mouse <- "8363"
    #  # expect_with_multi.mouse <- data.frame(
    #  #  ENTREZID.sp1 = "8363",
    #  #  ENTREZID.sp2 = c("102641229", "319155", "326620"),
    #  #  stringsAsFactors = FALSE
    #  #  )

    #  # Human gene that maps to multiple mouse ensembl gene ids via a single
    #  #   human ensembl gene id
    #  gene_ids_with_multi_mouse <- "26"
    #  expect_with_multi.mouse <- data.frame(
    #    ENTREZID.sp1 = rep("26", 4),
    #    ENTREZID.sp2 = c("243376", "243377", "69761", "76507"),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_homologues_with_biomart(
    #      gene_ids_with_multi_mouse
    #    ),
    #    expected = expect_with_multi.mouse,
    #    info = paste(
    #      "Human entrez id that maps to multiple mouse entrez ids",
    #      "(via a single HS-ensembl and multiple MM-ensembl ids)"
    #    )
    #  )

    #  # Human genes that map to the same mouse orthologues
    #  # "7314"/UBB and "7316"/UBC map to both "22187"/Ubb and "22190"/Ubc
    #  #
    #  # The test that uses UBB/UBC no longer works as of 14/12/2015
    #  #   since 7314 maps to a single gene and 7316 maps to a single gene now
    #  # gene_ids_many_to_many <- c("7314", "7316")
    #  # expect.many_to_many <- data.frame(
    #  #   ENTREZID.sp1 = rep(gene_ids_many_to_many, each = 2),
    #  #   ENTREZID.sp2 = rep(c("22187", "22190"), 2),
    #  #   stringsAsFactors = FALSE
    #  #   )
    #  # Hence, I rewrote the method using the human alkaline phosphatases ALPP
    #  #   and ALPPL2
    #  # Each of these genes maps to 3 different mouse gene ids (in biomart on
    #  #   14/12/15)
    #  gene_ids_many_to_many <- c("250", "251")
    #  expect.many_to_many <- data.frame(
    #    ENTREZID.sp1 = rep(c("250", "251"), each = 3),
    #    ENTREZID.sp2 = rep(c("11648", "11650", "76768"), times = 2),
    #    stringsAsFactors = FALSE
    #  )
    #  expect_equal(
    #    object = map_to_homologues_with_biomart(
    #      gene_ids_many_to_many
    #    ),
    #    expected = expect.many_to_many,
    #    info = "Genes that map to the same multiple genes"
    #  )
  }
)

###############################################################################

test_that(
  "map_to_homologues", {
    #
    #  # Map ensembl to ensembl: human to mouse, single one-to-one example
    #  expect_equal(
    #    object = map_to_homologues(
    #      gene_ids = "ENSG00000134294",
    #      sp1 = "hsapiens",
    #      sp2 = "mmusculus",
    #      idtype_sp1 = "ensembl_gene_id",
    #      idtype_sp2 = "ensembl_gene_id"
    #    ),
    #    expected = data.frame(
    #      ID.sp1 = "ENSG00000134294",
    #      ID.sp2 = "ENSMUSG00000022462",
    #      stringsAsFactors = FALSE
    #    ),
    #    info = paste("Ensembl-to-ensembl homologues; human to mouse singly",
    #                 "mapping example")
    #  )
    #  # Then reversed:
    #  expect_equal(
    #    object = map_to_homologues(
    #      gene_ids = "ENSMUSG00000022462",
    #      sp1 = "mmusculus",
    #      sp2 = "hsapiens",
    #      idtype_sp1 = "ensembl_gene_id",
    #      idtype_sp2 = "ensembl_gene_id"
    #    ),
    #    expected = data.frame(
    #      ID.sp1 = "ENSMUSG00000022462",
    #      ID.sp2 = "ENSG00000134294",
    #      stringsAsFactors = FALSE
    #    ),
    #    info = paste("Ensembl-to-ensembl homologues; mouse to human singly",
    #                 "mapping example")
    #  )
    #  # Human ensembl gene with no mouse orthologue
    #  # - note that default mapping is ensembl:ensembl and human:mouse
    #  # - ENSG00000284192 was found by searching in biomart with filters =
    #  # "with_mmusculus_homolog", values = FALSE and then keeping only those
    #  # genes with an entrez gene id
    #  expect_equal(
    #    object = map_to_homologues("ENSG00000284192"),
    #    expected = data.frame(
    #      ID.sp1 = "ENSG00000284192",
    #      ID.sp2 = as.character(NA),
    #      stringsAsFactors = FALSE
    #    ),
    #    info = "Non-mapping human gene, ensembl to ensembl"
    #  )
    #  # Human ensembl gene with multiple mouse orthologues
    #  expect_equal(
    #    object = map_to_homologues("ENSG00000002726"),
    #    expected = data.frame(
    #      ID.sp1 = "ENSG00000002726",
    #      ID.sp2 = c(
    #        "ENSMUSG00000029811",
    #        "ENSMUSG00000029813",
    #        "ENSMUSG00000039215",
    #        "ENSMUSG00000068536"
    #      ),
    #      stringsAsFactors = FALSE
    #    ),
    #    info = "One-to-many mapping human gene, ensembl-to-ensembl"
    #  )

    #  # Would like a gene that maps many:many for use in testing one.to.one =
    #  # TRUE
  }
)

###############################################################################

test_that(
  "map_to_homologues: one_to_one", {
    #  # One-to-one ensembl to ensembl: human to mouse, single one-to-one
    #  # example
    #  expect_equal(
    #    object = map_to_homologues(
    #      gene_ids = "ENSG00000134294",
    #      one_to_one = TRUE
    #    ),
    #    expected = data.frame(
    #      ID.sp1 = "ENSG00000134294",
    #      ID.sp2 = "ENSMUSG00000022462",
    #      stringsAsFactors = FALSE
    #    ),
    #    info = paste(
    #      "one_to_one ensembl-to-ensembl homologues; singly mapping example"
    #    )
    #  )
    #  expect_equal(
    #    object = map_to_homologues(
    #      gene_ids = "ENSG00000002726",
    #      one_to_one = TRUE
    #    ),
    #    expected = DF(
    #      ID.sp1 = "ENSG00000002726",
    #      ID.sp2 = as.character(NA)
    #    ),
    #    info = "one_to_one ensembl-to-ensembl homologues; 1:many example"
    #  )
    #  expect_equal(
    #    # Many to one example; note *29811 is one of four mouse homologues of
    #    # ENSG0000002726
    #    object = map_to_homologues(
    #      gene_ids = "ENSMUSG00000029811",
    #      sp1 = "mmusculus",
    #      sp2 = "hsapiens",
    #      one_to_one = TRUE
    #    ),
    #    expected = DF(
    #      ID.sp1 = "ENSMUSG00000029811",
    #      ID.sp2 = as.character(NA)
    #    ),
    #    info = "one_to_one ensembl-to-ensembl homologues; many:1 example"
    #  )
  }
)

###############################################################################

test_that("get_ensembl_homologue_field", {
  expect_equal(
    get_ensembl_homologue_field("hsapiens"),
    "hsapiens_homolog_ensembl_gene",
    info = "correct formatting of homology field for a given species"
  )
})

###############################################################################

test_that("is_valid_mart", {
  testthat::skip_if_not(requireNamespace("mockery", quietly = TRUE))

  # mock object for use in tests
  my_mart <- data.frame(required_id = 1:3)
  attr(my_mart, "class") <- "Mart"

  testthat::with_mock(
    listAttributes = mockery::mock(data.frame("required_id", "description")),
    expect_true(
      is_valid_mart(my_mart, "required_id"),
      "requested column is present in a Mart"
    ),
    .env = "biomaRt"
  )
})

###############################################################################
