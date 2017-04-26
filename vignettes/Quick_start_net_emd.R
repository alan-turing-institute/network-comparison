## ---- fig.show='hold'----------------------------------------------------
library("netdist")
# Set source directory and file properties for Virus PPI graph edge files
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
edge_format = "ncol"
file_pattern = ".txt"

# Calculate graphlet-based degree distributions for all orbits in graphlets 
# comprising up to 4 nodes for all graphs. This only needs to be done once 
# per graph (feature_type = "orbit", max_graphlet_size = 4).. 
# If feature_type is set to "feature_type", orbit counts for orbits in the
# same graphlet will be summed to generate graphlet counts
# If max_graphlet_size is set to 5, graphlet-based degree distributions will  
# be calculated for graphlets comprising up to 5 nodes.
virus_gdds <- gdd_for_all_graphs(
  source_dir = source_dir, format = edge_format, pattern = file_pattern, 
  feature_type = "orbit", max_graphlet_size = 4)
names(virus_gdds)

# Compute NetEMDs between all virus PPI graphs based on the computed graphlet- 
# based degree distributions using the default fast "optimise" method and no
# smoothing (default). The "optimise" method uses the built-in R optimise
# function to efficiently find the offset with the minimum EMD, but is not
# guaranteed to find the global minimum if EMD as a function of offset
# is non-convex and/or multimodal. The smoothing window width determines 
# whether to calculate the NetEMD from the unaltered discrete GDD histograms
# (smoothing_window_width = 0; default) or to first apply "nearest neighbour" 
# smoothing by "smearing" the discrete GDD histogram point masses across bins 
# of unit width (smoothing_window_width = 1). Returns a named list containing:
# (i) the NetEMDs and (ii) a table containing the graph names and indices 
# within the input GDD list for each pair of graphs compared.
out <- net_emds_for_all_graphs(virus_gdds, smoothing_window_width = 0)
print(out$net_emds)

# You can also specify method = "fixed_step" to use the much slower method of 
# exhaustively evaluating the EMD at all offsets separated by a fixed step. 
# The default step size is 1/2 the the minimum spacing between locations in 
# either histogram after normalising to unit variance. However, you can 
# specifiy your own fixed step using the optional "step_size" parameter.
# Note that this step size is applied to the histograms after they have been 
# normalised to unit variance

# Display NetEMDs for all network pairs, alongside the details of the network pairs
print(cbind(out$comp_spec, out$net_emds))

# Confirm NetEMD of a graph with itself is zero within (or at least close to)
# the limits of machine precision
purrr::map2(virus_gdds, virus_gdds, net_emd)
.Machine$double.eps

# The gdd_for_all_graphs and net_emds_for_all_graphs functions will run in 
# parallel using multiple threads where supported. The number of threads
# used is determined by the global R option "mc.cores". You can inspect the 
# current value of this using options("mc.cores") and set it with 
# options("mc.cores" = <num_cores>). To fully utilise a modern consumer
# processor, this should be set to 2x the number of available processor 
# cores as each core supports two threads.

