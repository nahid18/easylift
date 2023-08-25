library(GenomicRanges)
library(IRanges)

test_that("easylift function tests with valid chain files", {

  # Test 1: Test with a valid chain gzipped file
  chain_path_gz <- system.file("extdata", "hg19ToHg38.over.chain.gz", package = "easylift")
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 200),
    strand = "+"
  )
  genome(gr) <- "hg19"
  expect_type(easylift(gr, to = "hg38", chain = chain_path_gz), "S4")
  expect_no_error(easylift(gr, to = "hg38", chain = chain_path_gz))

  # Test 2: Test with a valid chain file
  chain_path <- system.file("extdata", "hg19ToHg38.over.chain", package = "easylift")
  gr2 <- GenomicRanges::GRanges(
    seqnames = "chr2",
    ranges = IRanges::IRanges(start = 200, end = 300),
    strand = "+"
  )
  genome(gr2) <- "hg19"
  expect_type(easylift(gr2, to = "hg38", chain = chain_path), "S4")
  expect_no_error(easylift(gr2, to = "hg38", chain = chain_path))

})

test_that("easylift function tests with error cases", {

  # Test 3: Test with an empty GRanges object
  gr3 <- GenomicRanges::GRanges()
  expect_error(easylift(gr3, to = "hg38", chain = chain_path))

  # Test 4: Test with missing genome information
  gr4 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  expect_error(easylift(gr4, to = "hg38", chain = chain_path))

})
