#' fastq html summary report
#'
#' Generate html summary report for the fastq files in a folder
#'
#' @param inpath A path to folder containing fastq/fastq.gz files
#' @return
#' @export

fq_report <- function(inpath) {
  calldir <- getwd()
  template <- paste0(system.file(package = 'rfaster2') , '/report.Rmd')

  rmarkdown::render(
    input = template,
    output_file = 'rfaster2-report.html',
    output_dir = calldir,
    params = list(
      fastq_dir = inpath
    )
  )
}
