## ------------------------------------------------------------------------
library("netdist")
# Set source directory and file properties for Virus PPI graph edge files
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
edge_format = "ncol"
file_pattern = "*"

# Load all graphs in the source folder matching the filename pattern
graphs <- read_simple_graphs(source_dir, format = edge_format, 
                             pattern = file_pattern)
print(names(graphs))

## ------------------------------------------------------------------------
ref_graph <- head(graphs, 1)[[1]]
query_graphs <- tail(graphs,- 1)

## ------------------------------------------------------------------------
# Set the maximum graphlet size to compute counts for
max_graphlet_size <- 4
neighbourhood_size <- 2

## ------------------------------------------------------------------------
expected_count_fn <- netdis_expected_graphlet_counts_ego_fn(
  ref_graph, max_graphlet_size, neighbourhood_size) 

## ------------------------------------------------------------------------
centred_counts <- purrr::map(query_graphs, netdis_centred_graphlet_counts,
                             max_graphlet_size = max_graphlet_size, 
                             neighbourhood_size = neighbourhood_size,
                             expected_ego_count_fn = expected_count_fn)

## ------------------------------------------------------------------------
# Netdis measure for graphlets of size 3
res3 <- netdis_for_all_graphs(centred_counts, 3)
netdis3_mat <- cross_comp_to_matrix(res3$netdis, res3$comp_spec)
# Netdis measure for graphlets of size 4
res4 <- netdis_for_all_graphs(centred_counts, 3)
netdis4_mat <- cross_comp_to_matrix(res4$netdis, res4$comp_spec)

## ------------------------------------------------------------------------
graphdists<-as.dist(netdis4_mat)
par(mfrow=c(1,2))
cex=1
# Dendrogram based on Netdis measure for graphlets of size 3
title = paste("Netdis: graphlet size = ", 3, sep = "")
plot(phangorn::upgma(as.dist(netdis3_mat), method="average"), use.edge.length=FALSE, 
     edge.width=cex*2, main=title, cex.lab=cex, cex.axis=cex, cex.main=cex, 
     cex.sub=cex, cex=cex)
# Dendrogram based on Netdis measure for graphlets of size 4
title = paste("Netdis: graphlet size = ", 4, sep = "")
plot(phangorn::upgma(as.dist(netdis4_mat), method="average"), use.edge.length=FALSE, 
     edge.width=cex*2, main=title, cex.lab=cex, cex.axis=cex, cex.main=cex, 
     cex.sub=cex, cex=cex)


