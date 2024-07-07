# rfaster2

> Get summary metrics about fastq data in `R`.

This package exports some of the functions of the Rust program [`faster2`]() in `R`. It uses compiled Rust code (thanks to the [`extendr`]() library) and is therefore pretty fast.

## Install

## Use

All `rfaster2` functions have the `fq_` prefix. All of them take a path to a fastx file (fastq or fastq.gz) and return various metrics - reads, bases, N50, quality scores...
