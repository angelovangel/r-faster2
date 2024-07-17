spk_tool <- function(label, x, values) {
  htmlwidgets::JS(
    sprintf(
      "function(sparkline, options, field){ return %s[field[0].offset]; }",
      jsonlite::toJSON(paste(label, x, ":",values, sep = " "))
    )
  )
}

# format the tooltip numbers back to their values
log_formatter <- htmlwidgets::JS(sprintf("function(x){ return Math.round(Math.pow(10, x)); }"))


### LENGTHS ###
len_density <- function(file) {

  density_obj <- rfaster2::fq_lengths(file) %>%
    as.numeric() %>%
    log10() %>%
    density(from = 2, to = 5, n = 100, na.rm = TRUE)
  density_obj
}

spark_len <- function(len_density_obj) {
  sparkline::spk_chr(paste( len_density_obj$x, ":", len_density_obj$y, sep = ""),
          #type = "bar",
          lineWidth = 3,
          fillColor = "#D0D3D4",
          lineColor = "#5D6D7E",
          spotColor = FALSE,
          minSpotColor = FALSE,
          maxSpotColor = "red",
          spotRadius = 3,
          width = 260, height = 40,
          tooltipFormat = "<span style='color: {{color}}'>&#9679;</span> {{prefix}}length: {{x}} {{suffix}}</span>",
          numberFormatter = log_formatter
  )
}
### LENGTHS ###

### GC-content
gc_density <- function(infile, nth) {

  density_obj <- rfaster2::fq_gc(infile, nth) %>%
    as.numeric() %>%
    # actually use density() here, not hist(). It returns a density list object with x and y, x is fixed from 1 to 100
    density(from = 0, to = 1, n = 100, na.rm = TRUE) # n is the number of equally spaced points at which the density is to be estimated.
  # return
  density_obj
}

spark_gc <- function(gc_density_obj) {
  spk_chr(paste( round(gc_density_obj$x, digits = 2), ":", gc_density_obj$y, sep = "" ),
          lineWidth = 3,
          fillColor = "#D0D3D4",
          lineColor = "#5D6D7E",
          spotColor = FALSE,
          minSpotColor = FALSE,
          maxSpotColor = "red",
          spotRadius = 3,
          width = 120, height = 40,
          tooltipFormat = "<span style='color: {{color}}'>&#9679;</span> {{prefix}}avg GC: {{x}} {{suffix}}</span>"
  )
}
###

### Q-SCORES ###

qscore_density <- function(file, phred, nth) {

  qscores <- rfaster2::fq_quals(infile = file, phred = phred, nth = nth)
  # actually use density() here, not hist(). It returns a density list object with x and y, x is fixed from 1 to 50
  density_obj <- density(qscores, from = 1, to = 60, n = 60, na.rm = TRUE) # n is the number of equally spaced points at which the density is to be estimated.
  density_obj
}


spark_phred <- function(phred_density_obj) {
  #spk_chr(paste( round(phred_density_obj$x, digits = 2), ":", phred_density_obj$y, sep = ""),
  fillcolor <- "#5D6D7E"
  sparkline::spk_chr(round(phred_density_obj$y, digits = 2),
                     type = "bar",
                     # to highlight q-value of 30, only array (60 elements) seems to work, don't know how to pass range map here
                     colorMap = c(
                       rep(fillcolor, 9), "red",
                       rep(fillcolor, 9), "red",
                       rep(fillcolor, 9), "red",
                       rep(fillcolor, 9), "red",
                       rep(fillcolor, 19)
                     ),
                     width = 140, height = 40,
                     tooltipFormatter = spk_tool("qscore ",phred_density_obj$x, round(phred_density_obj$y, 2))
  )
}
### Q-SCORES ###

### K-MERS ###
get_kmers <- function(infile, kmer, nth) {
  rfaster2::fq_kmers(infile, kmer, nth)
}

spark_kmers <- function(kmers_vec) {
  fillcolor <- "#5D6D7E"
  spk_chr(kmers_vec, type = "bar",
          barColor = "#5D6D7E",
          #colorMap = c("red", rep(fillcolor, 15), "red", rep(fillcolor, 15), "red", rep(fillcolor, 15), "red", rep(fillcolor, 15)),
          width = 220, height = 40
          #tooltipFormatter = spk_tool("", kmers_tbl$kmer, kmers_tbl$counts)
  )
}

### K-MERS ###
