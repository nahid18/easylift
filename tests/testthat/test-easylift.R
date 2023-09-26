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

  # Test 2 and 3: Test with a valid chain file
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

  # Test 4: Test with an empty GRanges object
  gr3 <- GenomicRanges::GRanges()
  expect_error(easylift(gr3, to = "hg38", chain = chain_path))

  # Test 5: Test with missing genome information
  gr4 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  expect_error(easylift(gr4, to = "hg38", chain = chain_path))
})

test_that("easylift succeeds with BiocFileCache", {
  # Load package
  library("easylift")

  # Create a test chain file in the temporary directory
  chain_file <- system.file("extdata", "hg19ToHg38.over.chain.gz", package = "easylift")

  # Test 6: Check if the chain file exists
  expect_true(file.exists(chain_file), "Chain file should exist.")

  # Create a BiocFileCache instance (to a temp location for this test script)
  path <- tempfile()
  bfc <- BiocFileCache(path, ask = FALSE)

  # Add the test chain file to the cache
  bfcadd(bfc, chain_file)

  # Query the cache for the chain file
  q <- bfcquery(bfc, chain_file)

  # Test 7: Check if the query result has rows (file exists in the cache)
  expect_true(nrow(q) > 0, "Chain file should exist in cache.")

  # Create a test GRanges object
  gr <- GenomicRanges::GRanges(
    seqname = Rle(c("chr1", "chr2"), c(100000, 100000)),
    ranges = IRanges::IRanges(start = 1, end = 200000)
  )

  genome(gr) <- "hg19"

  # Perform easylift with the target assembly
  tryCatch({
    result <- easylift(x=gr, to="hg38", bfc=bfc)
  }, error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    stop("easylift encountered an error.")
  })

  # Test 8: Check if easylift succeeded without error
  expect_true(!is(result, "try-error"), "easylift should succeed without error.")
})
