netemd_speed_test <- function() {
  ## load the data
  source_dir <- system.file(file.path("extdata", "random"), package = "netdist")
  print(source_dir)
  edge_format <- "ncol"
  file_pattern <- ""
  graphs <- read_simple_graphs(
    source_dir = source_dir,
    format = edge_format,
    pattern = file_pattern
  )
  n1 <- names(graphs)
  lab1 <- c()
  gdd_build_time <- c()
  netemd_time <- c()
  for (i in 1:length(graphs))
  {
    for (j in 1:(i))
    {
      g1 <- graphs[[i]]
      g2 <- graphs[[j]]
      lab1 <- append(lab1, paste(n1[i], n1[j], sep = ","))
      print(paste(n1[i], n1[j], sep = ","))
      fulltime_start <- Sys.time()
      gdd1 <- gdd(g1)
      gdd2 <- gdd(g2)
      netemd_start <- Sys.time()
      netemd_single_pair(gdd1, gdd2)
      end_time <- Sys.time()
      gdd_build_time <- append(
        gdd_build_time,
        as.double(netemd_start - fulltime_start)
      )
      netemd_time <- append(netemd_time, as.double(end_time - netemd_start))
    }
  }
  list(gdd_build_time = gdd_build_time, netemd_time = netemd_time)
}

#' @export
net_emd_speed_test_smooth <- function() {
  ## load the data
  source_dir <- system.file(file.path("extdata", "random"), package = "netdist")
  print(source_dir)
  edge_format <- "ncol"
  file_pattern <- ""
  graphs <- read_simple_graphs(
    source_dir = source_dir,
    format = edge_format,
    pattern = file_pattern
  )
  n1 <- names(graphs)
  lab1 <- c()
  gdd_build_time <- c()
  netemd_time <- c()
  for (i in 1:length(graphs))
  {
    for (j in 1:(i))
    {
      g1 <- graphs[[i]]
      g2 <- graphs[[j]]
      lab1 <- append(lab1, paste(n1[i], n1[j], sep = ","))
      print(paste(n1[i], n1[j], sep = ","))
      fulltime_start <- Sys.time()
      gdd1 <- gdd(g1)
      gdd2 <- gdd(g2)
      netemd_start <- Sys.time()
      net_emd(gdd1, gdd2, smoothing_window_width = 1)
      endTime <- Sys.time()
      gdd_build_time <- append(
        gdd_build_time,
        as.double(netemd_start - fulltime_start)
      )
      netemd_time <- append(netemd_time, as.double(endTime - netemd_start))
    }
  }
  list(gdd_build_time = gdd_build_time, netemd_time = netemd_time)
}
