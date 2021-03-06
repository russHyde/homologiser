% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homologue_db.R
\name{get_ensembl_homologue_field}
\alias{get_ensembl_homologue_field}
\title{get_ensembl_homologue_field}
\usage{
get_ensembl_homologue_field(sp)
}
\arguments{
\item{sp}{A species-name contraction, as used in
ensembl / `biomaRt`, for example, mmusculus for mouse and hsapiens for
human.}
}
\description{
In an ensembl biomart object, the column that contains the ensembl-gene ids
  for homologues in species `sp` is <sp>_homolog_ensembl_gene, where <sp>
  is the interpolated value of `sp` (eg, mmusculus).
}
