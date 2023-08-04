library("GenomeInfoDb")
library("rtracklayer")
library("R.utils")
library("tools")

#' Lift genomic coordinates from one genome assembly to another.
#'
#' This function takes a GRanges object with genomic coordinates in one genome
#' assembly and lifts them to another genome assembly using a chain file.
#'
#' @param gr A GRanges object with genomic coordinates in the original assembly.
#' @param to The target genome assembly (e.g., "hg38").
#' @param chain The path to the chain file containing the liftover mapping.
#' Can be provided in gzipped or non-gzipped format.
#' @return A GRanges object with lifted genomic coordinates.
#' @export
easylift <- function(gr, to, chain) {
  if (is.na(GenomeInfoDb::genome(gr))) {
    stop(
      "Error: The genome information for the 'gr' object is missing (NA). Please set the genome before using easylift."
    )
  }
  # Convert the input GRanges to the "UCSC" seqlevels style if not already
  if (GenomeInfoDb::seqlevelsStyle(gr) != "UCSC") {
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  }

  # Check if the chain file is gzipped and unzip if needed
  if (tools::file_ext(chain) == "gz") {
    tmp_dir <- tempdir()
    unzipped_chain <- file.path(tmp_dir, "tmp.chain")
    R.utils::gunzip(chain, destname=unzipped_chain, overwrite=TRUE, remove=FALSE)
    chain <- unzipped_chain
  }

  # Load the chain file
  ch <- rtracklayer::import.chain(chain)

  # Check if the 'to' genome is valid and available
  to_seqinfo <- GenomeInfoDb::getChromInfoFromUCSC(to, assembled.molecules.only = TRUE, as.Seqinfo = TRUE)
  if (is.null(to_seqinfo)) {
    stop(
      paste("Error: The genome", to, "is not available or recognized.")
    )
  }

  # LiftOver the genomic coordinates
  cur <- unlist(rtracklayer::liftOver(gr, ch))

  GenomeInfoDb::seqlevels(cur) <- GenomeInfoDb::seqlevels(to_seqinfo)
  GenomeInfoDb::seqinfo(cur) <- to_seqinfo

  return(cur)
}
