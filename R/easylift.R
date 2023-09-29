#' Lift genomic coordinates from one genome assembly to another.
#'
#' This function takes a GRanges object with genomic coordinates in one genome
#' assembly and lifts them to target genome assembly using a chain file.
#'
#' @param x A GRanges object with genomic coordinates in the original assembly.
#' @param to The target genome assembly (e.g., "hg38").
#' @param chain The path to the chain file containing the liftover mapping.
#' Can be provided in gzipped or non-gzipped format. If omitted, the function
#' will look in the default BiocFileCache for a properly named chain file.
#' @param bfc A BiocFileCache object (optional), if not provided (most typically)
#' the default location will be used.
#'
#' @return A GRanges object with lifted genomic coordinates.
#'
#' @examples
#' # Lift over the coordinates of the first 10 genes in the hg19 assembly
#' # to the hg38 assembly
#' library(easylift)
#' gr <- GRanges(
#'   seqname = Rle(c("chr1", "chr2"), c(100000, 100000)),
#'   ranges = IRanges(start = 1, end = 200000)
#' )
#' # Here, "hg19" is the source genome
#' genome(gr) <- "hg19"
#'
#' # Here, we use the `system.file()` function because the chain file is in the
#' # package (however if you need to point to any other file on your machine,
#' # just do 'chain <- "path/to/your/hg19ToHg38.over.chain.gz"'):
#' chain <- system.file("extdata", "hg19ToHg38.over.chain.gz", package = "easylift")
#'
#' # Here, "hg38" is the target genome
#' easylift(gr, "hg38", chain)
#'
#' \donttest{
#' # To use `BiocFileCache` for the chain file, add it to the cache as follows:
#' chain_file <- "/path/to/your/hg19ToHg38.over.chain.gz"
#' bfc <- BiocFileCache()
#'
#' # Add chain file to cache if already not available
#' if (nrow(bfcquery(bfc, basename(chain_file))) == 0)
#'    bfcadd(bfc, chain_file)
#'
#' # Then, use it in `easylift` like this:
#'
#' easylift(gr, "hg38")
#' # or
#' gr |> easylift("hg38")
#' }
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import rtracklayer
#' @import R.utils
#' @import tools
#' @import BiocFileCache
#' @import methods
#' @seealso
#' \code{\link[rtracklayer]{liftOver}} function from the \code{rtracklayer} package,
#' which is the basis for \code{easylift}.
#' @export
easylift <- function(x, to, chain, bfc) {
  args <- .perform_sanity_check(x, to, chain, bfc)
  x <- args$x
  to <- args$to
  chain <- args$chain
  bfc <- args$bfc
  unique_genomes <- args$unique_genomes

  # Convert the input GRanges to the "UCSC" seqlevels style if not already
  if (GenomeInfoDb::seqlevelsStyle(x) != "UCSC") {
    GenomeInfoDb::seqlevelsStyle(x) <- "UCSC"
  }

  if (.is_faulty(chain)) {
    chain <- .get_chain_from_BiocFileCache(unique_genomes, to, bfc)
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

### 'from' and 'to' are single strings of UCSC genome names, and 'bfc' is a BiocFileCache object
.get_chain_from_BiocFileCache <- function(from, to, bfc) {
  is_bfc_faulty <- .is_faulty(bfc)
  if (is_bfc_faulty) {
    bfc <- BiocFileCache()
  }
  capTo <-
    paste0(toupper(substr(to, 1, 1)), substr(to, 2, nchar(to)))
  trychainfile <- paste0(from, "To", capTo, ".over.chain")
  q <- bfcquery(bfc, trychainfile)
  if (nrow(q) == 0) {
    bfc_str <- if (is_bfc_faulty)
      "default"
    else
      "provided"
    stop(
      "Chain file not specified and filename with ",
      trychainfile,
      " pattern not found in BiocFileCache ",
      bfc_str,
      " location."
    )
  }
  chain <- bfc[[q$rid[1]]]
  return(chain)
}

# chain is expected to be string path of a chain file
.is_faulty <- function(y) {
  return(missing(y) || is.null(y))
}

.is_available <- function(y) {
  return(!missing(y) && !is.null(y))
}

.perform_sanity_check <- function(x, to, chain, bfc) {
  if (missing(chain)) {
    chain <- NULL
  }
  if (missing(bfc)) {
    bfc <- NULL
  }
  if (.is_faulty(x)) {
    stop("Please provide a 'GRanges' object 'x' with genomic coordinates.")
  }
  if (!is(x, "GRanges")) {
    stop("Please provide a valid 'GRanges' object 'x'.")
  }
  if (.is_faulty(to)) {
    stop("Please provide the target genome assembly 'to'.")
  }
  if (!is.character(to)) {
    stop("Please provide the target genome assembly 'to' as a string.")
  }
  if (.is_available(chain)) {
    if (!is.character(chain)) {
      stop("Please provide the chain file path 'chain' as a string.")
    }
  }
  if (.is_available(bfc)) {
    if (!is(bfc, "BiocFileCache")) {
      stop("Please provide a 'BiocFileCache' object 'bfc'.")
    }
  }
  if (anyNA(GenomeInfoDb::genome(x))) {
    stop("The genome information is missing. Please set genome(x) before using easylift.")
  }
  unique_genomes <- unique(GenomeInfoDb::genome(x))
  if (length(unique_genomes) > 1) {
    stop(
      "The 'GRanges' object 'x' contains genomic coordinates from multiple genomes. ",
      "Please provide 'x' with coordinates from a single genome assembly."
    )
  }
  return(list(
    x = x,
    to = to,
    chain = chain,
    bfc = bfc,
    unique_genomes = unique_genomes
  ))
}
