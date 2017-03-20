

# Virus Protein-Protein Interaction (PII) graphs
data_dir <- file.path("inst", "extdata", "VRPINS")

load_virus_data <- function(filename) {
  read_orca_graph(file = file.path(data_dir, filename), format = "ncol")
}

virusppi <- list(EBV = load_virus_data("EBV-1.txt"), 
                 ECL = load_virus_data("ECL-1.txt"),
                 HSV = load_virus_data("HSV-1-1.txt"),
                 KSHV = load_virus_data("KSHV-1.txt"),
                 VZV = load_virus_data("VZV-1.txt")
                 )

devtools::use_data(virusppi, overwrite = TRUE)

