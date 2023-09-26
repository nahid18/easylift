library(GenomicRanges)
library(IRanges)

test_that("easylift function tests with valid chain files", {
  # Test 1 and 2: Test with a valid chain gzipped file
  chain_path_gz <-
    system.file("extdata", "hg19ToHg38.over.chain.gz", package = "easylift")
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 200),
    strand = "+"
  )
  genome(gr) <- "hg19"
  expect_type(easylift(x = gr, to = "hg38", chain = chain_path_gz), "S4")
  expect_no_error(easylift(x = gr, to = "hg38", chain = chain_path_gz))

  # Test 3 and 4: Test with a valid chain file
  chain_path <-
    system.file("extdata", "hg19ToHg38.over.chain", package = "easylift")
  gr2 <- GenomicRanges::GRanges(
    seqnames = "chr2",
    ranges = IRanges::IRanges(start = 200, end = 300),
    strand = "+"
  )
  genome(gr2) <- "hg19"
  expect_type(easylift(x = gr2, to = "hg38", chain = chain_path), "S4")
  expect_no_error(easylift(x = gr2, to = "hg38", chain = chain_path))
})

test_that("easylift function tests with error cases", {
  # Test 5: Test with an empty GRanges object
  gr3 <- GenomicRanges::GRanges()
  expect_error(easylift(x = gr3, to = "hg38", chain = chain_path))

  # Test 6: Test with missing genome information
  gr4 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  expect_error(easylift(x = gr4, to = "hg38", chain = chain_path))

  # Test 7: Test with missing genome information
  gr5 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  genome(gr5) <- "hg19"
  expect_error(easylift(x = gr5, to = NULL, chain = chain_path))

  # Test 8: Test with missing genome information
  gr6 <- GenomicRanges::GRanges(
    seqnames = "chr4",
    ranges = IRanges::IRanges(start = 400, end = 500),
    strand = "+"
  )
  genome(gr6) <- "hg19"
  expect_error(easylift(x = gr6, to = NULL))
})

test_that("easylift succeeds with BiocFileCache", {
  # Load package
  library("easylift")

  # Test 9: Test if the chain file exists
  chain_file <-
    system.file("extdata", "hg19ToHg38.over.chain.gz", package = "easylift")
  expect_true(file.exists(chain_file), "Chain file should exist.")

  # Create a BiocFileCache instance (to a temp location for this test script)
  path <- tempfile()
  bfc <- BiocFileCache(path, ask = FALSE)

  # Add the test chain file to the cache
  bfcadd(bfc, chain_file)

  # Query the cache for the chain file
  q <- bfcquery(bfc, chain_file)

  # Test 10: Test if the query result has rows (file exists in the cache)
  expect_true(nrow(q) > 0, "Chain file should exist in cache.")

  # Create a test GRanges object
  gr <- GenomicRanges::GRanges(seqname = Rle(c("chr1", "chr2"), c(100000, 100000)),
                               ranges = IRanges::IRanges(start = 1, end = 200000))

  genome(gr) <- "hg19"

  # Test 11: Test success when bfc is provided
  tryCatch({
    result <- easylift(x = gr, to = "hg38", bfc = bfc)
  }, error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    stop("easylift encountered an error.")
  })

  expect_true(!is(result, "try-error"),
              "easylift should succeed without error.")

  # Test 12: Test success when bfc is NULL
  tryCatch({
    result2 <- easylift(x = gr, to = "hg38", bfc = NULL)
  }, error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    stop("easylift encountered an error.")
  })

  expect_true(!is(result2, "try-error"),
              "easylift should succeed without error.")

  # Test 13: Test success when chain is NULL
  tryCatch({
    result3 <- easylift(x = gr, to = "hg38", chain = NULL)
  }, error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    stop("easylift encountered an error.")
  })

  expect_true(!is(result3, "try-error"),
              "easylift should succeed without error.")

  # Test 14: Test success when both bfc and chain are NULL
  tryCatch({
    result4 <- easylift(x = gr, to = "hg38", chain = NULL, bfc = NULL)
  }, error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    stop("easylift encountered an error.")
  })

  expect_true(!is(result4, "try-error"),
              "easylift should succeed without error.")

  # Test 15: Test success when provided bfc is not a BiocFileCache instance
  # note: easylift will utilize default BiocFileCache location in this case
  tryCatch({
    result5 <- easylift(x = gr, to = "hg38", bfc = "bfc")
  }, error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    stop("easylift encountered an error.")
  })

  expect_true(!is(result5, "try-error"),
              "easylift should succeed without error.")

  # Test 16: Test error when chain file not found
  expect_error(easylift(x = gr, to = "hg100"))

  # Test 17: Test error when invalid chain file provided
  expect_error(easylift(x = gr, to = "hg100", chain="/random/path/not-found.chain"))
})
