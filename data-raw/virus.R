

# Virus Protein-Protein Interaction (PII) graphs
data_dir <- file.path("inst", "extdata", "VRPINS")

load_virus_data <- function(filename) {
  read_simple_graph(file = file.path(data_dir, filename), format = "ncol")
}

virusppi <- list(
  EBV = load_virus_data("EBV.txt"),
  ECL = load_virus_data("ECL.txt"),
  `HSV-1` = load_virus_data("HSV-1.txt"),
  KSHV = load_virus_data("KSHV.txt"),
  VZV = load_virus_data("VZV.txt")
)

devtools::use_data(virusppi, overwrite = TRUE)
