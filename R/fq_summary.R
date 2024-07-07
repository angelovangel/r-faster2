#' FASTQ file summary
#'
#' Obtain summary statistics for the records of a fastq/fastq.gz file
#'
#' @param infile A fastq/fastq.gz file
#' @return A named list
#' @export

fq_summary <- function(infile) {
  d <- rust_summary(infile)
  list(
    file = basename(infile),
    reads = d[1],
    bases = d[2],
    n_bases = d[3],
    min_len = d[4],
    max_len = d[5],
    N50 = d[6],
    gc_content = d[7]/d[2],
    q20 = d[8]/d[2],
    q30 = d[9]/d[2]
  )
}
