library("rtracklayer")
library("GenomeInfoDb")
library("tools")
library("R.utils")

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

# gr <- GRanges(seqname = Rle(paste("chr", 1, sep = "")),
#               ranges = IRanges(start = 1, end = 200000))
# genome(gr) <- "hg19"
# to <- "hg38"
# chain <- "data/hg19ToHg38.over.chain.gz"
# lifted <- easylift(gr, to, chain)
# print(seqinfo(gr))
# print(seqinfo(lifted))
