# rfaster2

> Get summary metrics about fastq data in `R`.

This package exports some of the functions of the Rust program [`faster2`]() in `R`. It uses Rust code (thanks to the [`extendr`]() library) and is therefore pretty fast.

## Install
If you have the rust toolchain:
```
devtools::install_github('angelovangel/rfaster2')
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

