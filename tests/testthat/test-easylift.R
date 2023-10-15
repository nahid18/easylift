library(GenomicRanges)
library(BiocFileCache)
library(IRanges)

test_that("easylift function tests with valid chain files", {
  # Test with a valid chain gzipped file
  chain_path <- system.file("extdata", "hg19ToHg38.over.chain", package = "easylift")
  gr1 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 200),
    strand = "+"
  )
  genome(gr1) <- "hg19"
  expect_type(easylift(x = gr1, to = "hg38", chain = chain_path), "S4")
  expect_no_error(easylift(x = gr1, to = "hg38", chain = chain_path))
})

test_that("easylift function tests with error cases", {
  # Test with an empty GRanges object
  gr3 <- GenomicRanges::GRanges()
  expect_error(easylift(x = gr3, to = "hg38", chain = chain_path))

  # Test with missing genome information
  gr4 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  expect_error(easylift(x = gr4, to = "hg38", chain = chain_path))

  # Test with missing genome information
  gr5 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  genome(gr5) <- "hg19"
  expect_error(easylift(x = gr5, to = NULL, chain = chain_path))

  # Test with missing genome information
  gr6 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  genome(gr6) <- "hg19"
  expect_error(easylift(x = gr6, to = NULL))
})

test_that("easylift succeeds with BiocFileCache", {
  # Create a BiocFileCache instance (to a temp location for this test script)
  path <- tempfile()
  bfc <- BiocFileCache(path, ask = FALSE)

  chain_file <- system.file("extdata", "hg19ToHg38.over.chain", package = "easylift")

  # Add the test chain file to the cache
  bfcadd(bfc, chain_file)

  # Query the cache for the chain file
  q <- bfcquery(bfc, chain_file)

  # Test if the query result has rows (file exists in the cache)
  expect_true(nrow(q) > 0, "Chain file should exist in cache.")

  # Create a test GRanges object
  gr <- GenomicRanges::GRanges(seqname = Rle(c("chr1", "chr2"), c(100000, 100000)),
                               ranges = IRanges::IRanges(start = 1, end = 200000))

  genome(gr) <- "hg19"

  # Should pass since bfc is provided properly to easylift
  expect_no_error(easylift(gr, "hg38", NULL, bfc))
  expect_no_error(gr |> easylift("hg38", NULL, bfc))
  expect_no_error(easylift(x = gr, to = "hg38", bfc = bfc))

  # Test error when provided bfc is not a BiocFileCache instance
  expect_error(easylift(gr, "hg38", bfc))
  expect_error(gr |> easylift("hg38", bfc))

  # Should error because the chain file was added to tempfile() but not in the
  # default cache location
  expect_error(easylift(gr, "hg38"))
  expect_error(gr |> easylift("hg38"))

  # Test error when provided bfc is not a BiocFileCache instance
  expect_error(easylift(x = gr, to = "hg38", bfc = "bfc"))

  # Test error when chain file not found
  expect_error(easylift(x = gr, to = "hg100"))

  # Test error when invalid chain file provided
  expect_error(easylift(x = gr, to = "hg100", chain="/random/path/not-found.chain"))
})
