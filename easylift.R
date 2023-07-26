library(rtracklayer)
library(GenomicRanges)

easylift <- function(gr, to, chain) {
  if (is.na(genome(gr))) {
    stop(
      "Error: The genome information for the 'gr' object is missing (NA). Please set the genome before using easylift."
    )
  }
  ch <- rtracklayer::import.chain(chain)
  seqlevelsStyle(gr) <- "UCSC"
  
  cur <- unlist(rtracklayer::liftOver(gr, ch))
  genome(cur) <- to
  
  to_seqinfo <-
    getChromInfoFromUCSC(to,
                         assembled.molecules.only = TRUE,
                         as.Seqinfo = TRUE)
  lifted_seqlevels <- seqlevels(cur)[genome(cur) %in% to]
  seqinfo(cur) <- update(seqinfo(cur), to_seqinfo[lifted_seqlevels])
  return(cur)
}

gr <- GRanges(seqname = Rle(paste("chr", 1, sep = "")),
              ranges = IRanges(start = 266038, end = 16207633))
genome(gr) <- "hg19"
to <- "hg38"
chain <- "data/hg19ToHg38.over.chain"

lifted <- easylift(gr, to, chain)

print(gr)
print(lifted)
print(seqlengths(lifted))
print(mcols(cur))
print(length(gr) - length(lifted))