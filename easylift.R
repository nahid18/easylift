library(rtracklayer)
library(GenomicRanges)

easylift <- function(gr, to, chain) {
  if (is.na(genome(gr))) {
    stop(
      "Error: The genome information for the 'gr' object is missing (NA). Please set the genome before using easylift."
    )
  }
  # Convert the input GRanges to the "UCSC" seqlevels style if not already
  if (seqlevelsStyle(gr) != "UCSC") {
    seqlevelsStyle(gr) <- "UCSC"
  }


  # Check if the chain file is gzipped and unzip if needed
  if (tools::file_ext(chain) == "gz") {
    tmp_dir <- tempdir()
    unzipped_chain <- file.path(tmp_dir, "unzipped.chain")
    gunzip(chain, destname=unzipped_chain, overwrite=TRUE, remove=FALSE)
    chain <- unzipped_chain
  }
  
  # Load the chain file
  ch <- import.chain(chain)
  
  # Check if the 'to' genome is valid and available
  to_genome_info <- getChromInfoFromUCSC(to, assembled.molecules.only = TRUE, as.Seqinfo = TRUE)
  if (is.null(to_genome_info)) {
    stop(
      paste("Error: The genome", to, "is not available or recognized.")
    )
  }
  
  # LiftOver the genomic coordinates
  cur <- unlist(liftOver(gr, ch))
  genome(cur) <- to
  
  # Update the seqinfo with the target genome information
  lifted_seqlevels <- seqlevels(cur)[genome(cur) %in% to]
  if (length(lifted_seqlevels) == 0) {
    stop("Error: LiftOver did not produce any valid results for the target genome.")
  }
  to_seqinfo <- to_genome_info[lifted_seqlevels]
  seqinfo(cur) <- update(seqinfo(cur), to_seqinfo)
  
  return(cur)
}

gr <- GRanges(seqname = Rle(paste("chr", 1, sep = "")),
              ranges = IRanges(start = 1, end = 200000))
genome(gr) <- "hg19"
to <- "hg38"
chain <- "data/hg19ToHg38.over.chain.gz"

lifted <- easylift(gr, to, chain)

print(gr)
print(lifted)
print(seqlengths(lifted))
print(mcols(lifted))
print(length(gr) - length(lifted))