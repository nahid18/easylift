#' Lift genomic coordinates from one genome assembly to another.
#'
#' This function takes a GRanges object with genomic coordinates in one genome
#' assembly and lifts them to another genome assembly using a chain file.
#'
#' @param x A GRanges object with genomic coordinates in the original assembly.
#' @param to The target genome assembly (e.g., "hg38").
#' @param chain The path to the chain file containing the liftover mapping.
#' Can be provided in gzipped or non-gzipped format.
#' @return A GRanges object with lifted genomic coordinates.
#' @examples
#' # Lift over the coordinates of the first 10 genes in the hg19 assembly
#' # to the hg38 assembly
#' library(GenomicRanges)
#' library(easylift)
#' gr <- GRanges(seqname = Rle(c("chr1", "chr2"), c(100000, 100000)), 
#' ranges = IRanges(start = 1, end = 200000))
#' genome(gr) <- "hg19"
#' chain <- system.file("extdata/hg19ToHg38.over.chain.gz", package="easylift")
#' easylift(gr, "hg38", chain)
#' @import GenomeInfoDb
#' @import rtracklayer
#' @import R.utils
#' @import tools
#' @import BiocFileCache
#' @export
easylift <- function(x, to, chain) {
  if (any(is.na(GenomeInfoDb::genome(x)))) {
    stop(
      "The genome information is missing. Please set genome(x) before using easylift."
    )
  }
  # Convert the input GRanges to the "UCSC" seqlevels style if not already
  if (GenomeInfoDb::seqlevelsStyle(x) != "UCSC") {
    GenomeInfoDb::seqlevelsStyle(x) <- "UCSC"
  }

  if (missing(chain)) {
    bfc <- BiocFileCache()
    capTo <- paste0(toupper(substr(to,1,1)),substr(to,2,nchar(to))) # capitalize first letter
    trychainfile <- paste0(genome(x),"To",capTo,".over.chain")
    q <- bfcquery(bfc, trychainfile)
    if (nrow(q) >= 1) {
      chain <- bfc[[q$rid[1]]]
    } else {
      stop("Chain file not specified and not found in BiocFileCache")
    }
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
      paste("The genome", to, "is not available or recognized.")
    )
  }

  # LiftOver the genomic coordinates
  cur <- unlist(rtracklayer::liftOver(x, ch))

  GenomeInfoDb::seqlevels(cur) <- GenomeInfoDb::seqlevels(to_seqinfo)
  GenomeInfoDb::seqinfo(cur) <- to_seqinfo

  return(cur)
}
