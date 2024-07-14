#' fastq html summary report
#'
#' Generate html summary report for the fastq files in a folder
#'
#' @param inpath A path to folder containing fastq/fastq.gz files
#' @param pattern Regex pattern to search for fastq files
#' @param platform Sequencing platform used
#' @param outfile name of html report file, default is rfaster2-report.html, written in the calling directory
#' @export

fq_report <- function(inpath, pattern = 'fast(q|q.gz)$', platform = 'Nanopore', outfile = 'rfaster2-report.html') {
  calldir <- getwd()
  template <- paste0(system.file(package = 'rfaster2') , '/report.Rmd')

  rmarkdown::render(
    input = template,
    output_file = outfile,
    output_dir = calldir,
    params = list(
      fastq_dir = inpath,
      fastq_pattern = pattern,
      platform = platform
    )
  )
}
