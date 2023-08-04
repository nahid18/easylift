
<!-- README.md is generated from README.Rmd. Please edit that file -->

# easylift

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/nahid18/easylift)](https://github.com/nahid18/easylift/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/nahid18/easylift)](https://github.com/nahid18/easylift/pulls)
<!-- badges: end -->

The goal of `easylift` is to perform genomic liftover given `GRanges`
and `chain` file.

## Installation instructions

<!-- Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `easylift` from [Bioconductor](http://bioconductor.org/) using the following code: -->
<!-- ```{r 'install', eval = FALSE} -->
<!-- if (!requireNamespace("BiocManager", quietly = TRUE)) { -->
<!--     install.packages("BiocManager") -->
<!-- } -->
<!-- BiocManager::install("easylift") -->
<!-- ``` -->

Get the latest development version from
[GitHub](https://github.com/nahid18/easylift) with:

``` r
BiocManager::install("nahid18/easylift")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("easylift")
library("GenomicRanges")
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
gr <- GRanges(seqname = Rle(paste("chr", 1, sep = "")),
              ranges = IRanges(start = 1, end = 200000))
genome(gr) <- "hg19"
to <- "hg38"
chain <- "hg19ToHg38.over.chain.gz"
lifted <- easylift(gr, to, chain)
```

## Citation

Below is the citation output from using `citation('easylift')` in R.
Please run this yourself to check for any updates on how to cite
**easylift**.

``` r
print(citation('easylift'), bibtex = TRUE)
```

Please note that the `easylift` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `easylift` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://nahid18.github.io/easylift) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
