netEMDSpeedTest <- function() {
  ## load the data
  source_dir <- system.file(file.path("extdata", "random"), package = "netdist")
  print(source_dir)
  edge_format <- "ncol"
  file_pattern <- ""
  #    source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
  #    edge_format = "ncol"
  #    file_pattern = ".txt"
  graphs <- read_simple_graphs(source_dir = source_dir, format = edge_format, pattern = file_pattern)
  n1 <- names(graphs)
  lab1 <- c()
  gddBuildTime <- c()
  netEMDtime <- c()
  for (i in 1:length(graphs))
  {
    for (j in 1:(i))
    {
      g1 <- graphs[[i]]
      g2 <- graphs[[j]]
      lab1 <- append(lab1, paste(n1[i], n1[j], sep = ","))
      print(paste(n1[i], n1[j], sep = ","))
      fulltimeStart <- Sys.time()
      gdd1 <- gdd(g1)
      gdd2 <- gdd(g2)
      netEMDStart <- Sys.time()
      net_emd(gdd1, gdd2)
      endTime <- Sys.time()
      gddBuildTime <- append(gddBuildTime, as.double(netEMDStart - fulltimeStart))
      netEMDtime <- append(netEMDtime, as.double(endTime - netEMDStart))
    }
  }
  list(gddBuildTime = gddBuildTime, netEMDtime = netEMDtime)
}
