# rfaster2

> Get summary metrics about fastq data in `R`.

This package exports some of the functions of the Rust program [`faster2`]() in `R`. It uses Rust code (thanks to the [`extendr`]() library) and is therefore pretty fast.

## Install

If you have the rust toolchain

``` r
remotes::install_github('angelovangel/rfaster2')
```

If you don't have it, Rust can be installed with

``` bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

If you don't want to setup Rust, use compiled GitHub releases. Binary packages are available for Linux (x86_64)

``` r
install.packages(
https://github.com/angelovangel/rfaster2/releases/download/v0.1.0/rfaster2_0.1.0_R_x86_64-pc-linux-gnu.tar.gz
repos = NULL
)
```

and for macOS (aarch64, a.k.a. Apple Silicon)

``` r
install.packages(
  "https://github.com/angelovangel/rfaster2/releases/download/v0.1.0/rfaster2_0.1.0_macos_aarch64.tgz",
  repos = NULL
)
```

## Use

All `rfaster2` functions have the `fq_` prefix. All of them take a path to a fastx file (fastq or fastq.gz) and return various metrics - reads, bases, N50, quality scores...

``` r
# get 'mean' qscores for the reads in a file
fq_quals(infile = path/to/fastqfile, phred = TRUE)

# get summary for a file
fq_summary(infile = path/to/fastqfile)

# get a data frame with statistics for a list of fastq files 
f <- list.files('path', pattern = 'fastq', full.names = T)
lapply(f, fq_summary) %>% dplyr::bind_rows()
```
