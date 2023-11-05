
# easylift

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/nahid18/easylift)](https://github.com/nahid18/easylift/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/nahid18/easylift)](https://github.com/nahid18/easylift/pulls)
<!-- badges: end -->

The goal of `easylift` is to perform genomic liftover given `GRanges`
and `chain` file. `easylift` was developed as part of the 
[tidyomics project](https://www.biorxiv.org/content/10.1101/2023.09.10.557072v2)'s
[open challenges](https://github.com/orgs/tidyomics/projects/1).

## Installation

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("easylift")
```

## Documentation

To view documentation:

``` r
browseVignettes("easylift")
```

## Import

Import the libraries

``` r
library("easylift")
```

Call `easylift` with `GRanges` object, target genome and the chain file.

``` r
gr <- GRanges(
  seqname = Rle(
    c("chr1", "chr2"), 
    c(100000, 100000)
  ),
  ranges = IRanges(
    start = 1, 
    end = 200000
  )
)
# Here, "hg19" is the source genome
genome(gr) <- "hg19"
chain <- "hg19ToHg38.over.chain.gz"

# Here, "hg38" is the target genome
easylift(gr, "hg38", chain)
```

### BiocFileCache

To use `BiocFileCache` for the chain file, add it to the cache:

``` r
chain_file <- "/path/to/your/hg19ToHg38.over.chain.gz"
bfc <- BiocFileCache()

# Add chain file to cache if already not available
if (nrow(bfcquery(bfc, basename(chain_file))) == 0)
    bfcadd(bfc, chain_file)
```

Then, use it in `easylift`:

``` r
easylift(gr, "hg38") 
# or
gr |> easylift("hg38") 
```

## Citation

To cite package `easylift` in publications use:

Al Nahid A, Pagès H, Love M (2023). *easylift: An R package to perform
genomic liftover*. R package version 1.0.0,
<https://github.com/nahid18/easylift>.

A BibTeX entry for LaTeX users is

      @Manual{,
        title = {easylift: An R package to perform genomic liftover},
        author = {Abdullah Al Nahid, Hervé Pagès, Michael Love},
        year = {2023},
        note = {R package version 1.0.0},
        url = {https://github.com/nahid18/easylift},
      }

Please note that the `easylift` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `easylift` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
