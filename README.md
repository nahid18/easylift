
# easylift

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/nahid18/easylift)](https://github.com/nahid18/easylift/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/nahid18/easylift)](https://github.com/nahid18/easylift/pulls)
<!-- badges: end -->

The goal of `easylift` is to perform genomic liftover given `GRanges`
and `chain` file.

## Installation

Get the latest development version from
[GitHub](https://github.com/nahid18/easylift) with:

``` r
BiocManager::install("nahid18/easylift")
```

## Usage

Import the libraries

``` r
library("easylift")
library("GenomicRanges")
```

Call `easylift` with `GRanges` object, target genome and chain file.

``` r
gr <- GRanges(seqname = Rle(paste("chr", 1, sep = "")),
              ranges = IRanges(start = 1, end = 200000))
genome(gr) <- "hg19"
to <- "hg38"
chain <- "hg19ToHg38.over.chain.gz"

easylift(gr, to, chain)
```

## Citation

To cite package ‘easylift’ in publications use:

Al Nahid A, Love M (2023). *easylift: An R package to perform genomic
liftOver*. R package version 0.0.3,
<https://github.com/nahid18/easylift>.

A BibTeX entry for LaTeX users is

      @Manual{,
        title = {easylift: An R package to perform genomic liftOver},
        author = {Abdullah Al Nahid and Michael Love},
        year = {2023},
        note = {R package version 0.0.3},
        url = {https://github.com/nahid18/easylift},
      }

Please note that the `easylift` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `easylift` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
