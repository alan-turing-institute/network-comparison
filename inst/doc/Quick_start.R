## ---- fig.show='hold'----------------------------------------------------
library("netdist")
# Load example virus PPI graphs from files
virus_edges <- read_all_graphs_as_orca_edge_lists(
  system.file(package = "netdist", "extdata", "VRPINS"),
  format = "ncol", pattern = ".txt")
attr(virus_edges, "names")

# Calculate graphlet orbit degree distributions up to 4 nodes for all graphs 
# This only needs to be done once per graph
virus_godd <- purrr::map(virus_edges, godd)

# Generate a cross-comparison matrix listing all combinations of graphs
comp_spec <- graph_cross_comparison_spec(virus_edges)
comp_spec[1:5,]

# Compute NetEMD between all virus PPI graphs based on all graphlet orbir
# degree distributions up to 4 nodes
net_emds <- purrr::simplify(
  purrr::map2(comp_spec$index_a, comp_spec$index_b, function(index_a, index_b) {
  net_emd(virus_godd[[index_a]], virus_godd[[index_b]])
}))

# Link NetEMDs with their respective comp_specs
comp_spec$net_emd = net_emds
print(comp_spec)

# Confirm NetEMD of a graph with itself is zero within (or at least close to)
# the limits of machine precision
purrr::map2(virus_godd, virus_godd, net_emd)
.Machine$double.eps

