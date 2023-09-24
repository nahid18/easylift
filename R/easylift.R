#' Lift genomic coordinates from one genome assembly to another.
#'
#' This function takes a GRanges object with genomic coordinates in one genome
#' assembly and lifts them to target genome assembly using a chain file.
#'
#' @param x A GRanges object with genomic coordinates in the original assembly.
#' @param to The target genome assembly (e.g., "hg38").
#' @param chain The path to the chain file containing the liftover mapping.
#' Can be provided in gzipped or non-gzipped format.
#' @return A GRanges object with lifted genomic coordinates.
#' @examples
#' # Lift over the coordinates of the first 10 genes in the hg19 assembly
#' # to the hg38 assembly
#' library(easylift)
#' gr <- GRanges(
#'   seqname = Rle(c("chr1", "chr2"), c(100000, 100000)),
#'   ranges = IRanges(start = 1, end = 200000)
#' )
#' genome(gr) <- "hg19"
#'
#' # here, we use the `system.file()` function because the chain file is in the
#' # package (however if you need to point to any other file on your machine,
#' # just do 'chain <- "path/to/your/file"'):
#' chain <- system.file("extdata/hg19ToHg38.over.chain.gz", package="easylift")
#' easylift(gr, "hg38", chain)
#' 
#' # To use `BiocFileCache` for the chain file, you can add it to the cache as follows:
#' # chain_file <- "/path/to/your/chainfile.chain"
#' # bfc <- BiocFileCache()
#' # bfcadd(bfc, chain_file)
#'
#' # Then, you can use it in `easylift` like this:
#' # easylift(gr, "hg38")
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import rtracklayer
#' @import R.utils
#' @import tools
#' @import BiocFileCache
#' @seealso
#' \code{\link[rtracklayer]{liftOver}} function from the \code{rtracklayer} package,
#' which is the basis for \code{easylift}.
#' @export
easylift <- function(x, to, chain) {
  # Check if genome(x) has multiple unique values
  unique_genomes <- unique(GenomeInfoDb::genome(x))
  if (length(unique_genomes) > 1) {
    stop("The 'GRanges' object 'x' contains genomic coordinates from multiple genomes. ",
         "Please provide 'x' with coordinates from a single genome assembly.")
  }

  if (anyNA(GenomeInfoDb::genome(x))) {
    stop("The genome information is missing. Please set genome(x) before using easylift.")
  }
  # Convert the input GRanges to the "UCSC" seqlevels style if not already
  if (GenomeInfoDb::seqlevelsStyle(x) != "UCSC") {
    GenomeInfoDb::seqlevelsStyle(x) <- "UCSC"
  }

  if (missing(chain)) {
    bfc <- BiocFileCache()
    capTo <-
      paste0(toupper(substr(to, 1, 1)), substr(to, 2, nchar(to))) # capitalize first letter
    trychainfile <- paste0(genome(x), "To", capTo, ".over.chain")
    q <- bfcquery(bfc, trychainfile)
    if (nrow(q) >= 1) {
      chain <- bfc[[q$rid[1]]]
    } else {
      stop(
        "Chain file not specified and not found in BiocFileCache. Look for a chain file named '",
        trychainfile,
        "' for liftover to '",
        to,
        "' in your BiocFileCache or provide the 'chain' argument with the file path."
      )
    }
  }

  # Check if the chain file is gzipped and unzip if needed
  if (tools::file_ext(chain) == "gz") {
    tmp_dir <- tempdir()
    unzipped_chain <- file.path(tmp_dir, "tmp.chain")
    R.utils::gunzip(chain,
                    destname = unzipped_chain,
                    overwrite = TRUE,
                    remove = FALSE)
    chain <- unzipped_chain
  }

  # Load the chain file
  ch <- rtracklayer::import.chain(chain)

  # Check if the 'to' genome is valid and available
  to_seqinfo <-
    GenomeInfoDb::getChromInfoFromUCSC(to,
                                       assembled.molecules.only = TRUE,
                                       as.Seqinfo = TRUE)
  if (is.null(to_seqinfo)) {
    stop("The genome", to, "is not available or recognized.")
  }

  # LiftOver the genomic coordinates
  cur <- unlist(rtracklayer::liftOver(x, ch))

  GenomeInfoDb::seqlevels(cur) <-
    GenomeInfoDb::seqlevels(to_seqinfo)
  GenomeInfoDb::seqinfo(cur) <- to_seqinfo

  return(cur)
}
