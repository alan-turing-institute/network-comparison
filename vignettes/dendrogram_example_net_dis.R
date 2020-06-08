## ------------------------------------------------------------------------
#A CHANGE TO FILE- COMMENT.
library("netdist")
edge_format = "ncol"
# Load reference graph (used for Netdis. Not required for NetEMD)
ref_path = file.path(system.file(file.path("extdata", "random"),
                                 package = "netdist"),
                     "ER_1250_10_1")
ref_graph <- read_simple_graph(ref_path, format = edge_format)

# Set source directory and file properties for Virus PPI graph edge files
source_dir <- system.file(file.path("extdata", "VRPINS"),
                          package = "netdist")
edge_format <- "ncol"
file_pattern <- "*"

# Load all graphs in the source folder matching the filename pattern
query_graphs <- read_simple_graphs(source_dir,
                                   format = edge_format, 
                                   pattern = file_pattern)
print(names(query_graphs))

## ------------------------------------------------------------------------
# Set the maximum graphlet size to compute counts for
max_graphlet_size <- 4
neighbourhood_size <- 2

## ------------------------------------------------------------------------

# Calculate netdis measure for graphlets up to size max_graphlet_size
netdis_result <- netdis_many_to_many(query_graphs,
                                     ref_graph,
                                     max_graphlet_size = max_graphlet_size,
                                     neighbourhood_size = neighbourhood_size)

# Netdis measure for graphlets of size 3
res3 <- netdis_result$netdis["netdis3", ]
netdis3_mat <- cross_comp_to_matrix(res3, netdis_result$comp_spec)

print("Netdis: graphlet size = 3")
print(netdis3_mat)

# Netdis measure for graphlets of size 4
res4 <- netdis_result$netdis["netdis4", ]
netdis4_mat <- cross_comp_to_matrix(res4, netdis_result$comp_spec)

print("Netdis: graphlet size = 4")
print(netdis4_mat)

## ------------------------------------------------------------------------
graphdists <- as.dist(netdis4_mat)
par(mfrow = c(1, 2))
cex <- 1

# Dendrogram based on Netdis measure for graphlets of size 3
title <- paste("Netdis: graphlet size = ", 3, sep = "")
plot(phangorn::upgma(as.dist(netdis3_mat), method = "average"),
     use.edge.length = FALSE, 
     edge.width = cex*2,
     main = title,
     cex.lab = cex, cex.axis = cex,
     cex.main = cex, cex.sub = cex,
     cex = cex)

# Dendrogram based on Netdis measure for graphlets of size 4
title = paste("Netdis: graphlet size = ", 4, sep = "")
plot(phangorn::upgma(as.dist(netdis4_mat), method = "average"),
     use.edge.length = FALSE, 
     edge.width = cex*2,
     main = title,
     cex.lab = cex, cex.axis = cex,
     cex.main = cex, cex.sub = cex,
     cex = cex)

## ------------------------------------------------------------------------
cex <- 1.5
col <- colorRampPalette(colors = c("blue","white"))(100)
title <- paste("Netdis: graphlet size = ", 3, sep = "")
heatmap(netdis3_mat, Rowv = NULL, Colv = NULL, col = col, main = title,
        cexRow = cex, cexCol = cex, symm = TRUE)

## ------------------------------------------------------------------------
cex <- 1.5
col <- colorRampPalette(colors = c("blue","white"))(100)
title <- paste("Netdis: graphlet size = ", 4, sep = "")
heatmap(netdis4_mat, Rowv = NULL, Colv = NULL, col = col, main = title,
        cexRow = cex, cexCol = cex, symm = TRUE)

