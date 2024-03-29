## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ---- packages, message= FALSE------------------------------------------------
# Load packages/libraries
library("netdist")
library("igraph")

## ---- graphs,fig.align='center'-----------------------------------------------
# Set source directory for Virus protein-protein interaction edge files stored in the netdist package.
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")


# Load query graphs as igraph objects
graph_1 <- read_simple_graph(file.path(source_dir, "EBV.txt"),
                             format = "ncol")

graph_2 <- read_simple_graph(file.path(source_dir, "ECL.txt"),
                             format = "ncol")

# Herpes virus EBV protein-protein interaction graph with 60 nodes and 208 edges.
graph_1

# Herpes virus ECL protein-protein interaction graph with 1941 nodes and 3989 edges.
graph_2

#A simple visualization of the graphs.
plot(graph_1,vertex.size=0.5,vertex.label=NA)
plot(graph_2,vertex.size=0.5,vertex.label=NA)

## ---- netemd,fig.align='center'-----------------------------------------------
# Set source directory for Virus protein-protein interaction network edge files stored in the netdist package.
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

# Load query graphs as igraph objects
# Herpes virus EBV protein-protein interaction graph with 60 nodes and 208 edges.
graph_1 <- read_simple_graph(file.path(source_dir, "EBV.txt"),
                             format = "ncol")

# Herpes virus ECL protein-protein interaction graph with 1941 nodes and 3989 edges.
graph_2 <- read_simple_graph(file.path(source_dir, "ECL.txt"),
                             format = "ncol")

# One to one NetEmd comparison.
netemd_one_to_one(graph_1=graph_1,graph_2=graph_2,feature_type="orbit",smoothing_window_width = 1)#Use of smoothing window 1 is given for discrete integer distributions. If the network features are considered continuous variables smoothing_window_width equal to zero is recommended.

## ---- netemdEigen,fig.align='center'------------------------------------------
#Laplacian
Lapg_1 <- igraph::laplacian_matrix(graph = graph_1,normalized = FALSE,sparse = FALSE)
Lapg_2 <- igraph::laplacian_matrix(graph = graph_2,normalized = FALSE,sparse = FALSE)

#Normalized Laplacian
NLapg_1 <- igraph::laplacian_matrix(graph = graph_1,normalized = TRUE,sparse = FALSE)
NLapg_2 <- igraph::laplacian_matrix(graph = graph_2,normalized = TRUE,sparse = FALSE)

#Spectra (This may take a couple of minutes).
props_1 <- cbind(L.Spectra= eigen(Lapg_1)$values, NL.Spectra= eigen(NLapg_1)$values) 
props_2 <- cbind(L.Spectra= eigen(Lapg_2)$values, NL.Spectra= eigen(NLapg_2)$values) 

netemd_one_to_one(dhists_1 = props_1,dhists_2 = props_2,smoothing_window_width = 0)#Use of smoothing window 1 is given for discrete integer distributions. If the network features are considered continuous variables smoothing_window_width equal to zero is recommended.

## ----netdisgoldstand,fig.align='center'---------------------------------------
# Lattice graphs to be used as gold-standard as a reference point comparison
goldstd_1 <- igraph::graph.lattice(c(8,8)) #Graph with 8^2 nodes
goldstd_2 <- igraph::graph.lattice(c(44,44)) #Graph with 44^2 nodes

plot(goldstd_1,vertex.size=0.8,vertex.label=NA)
plot(goldstd_2,vertex.size=0.5,vertex.label=NA)


# Netdis using the goldstd_1 graph as gold-standard reference point
netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = goldstd_1)

# Netdis using the goldstd_2 graph as gold-standard reference point
netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = goldstd_2)

## ---- netdisGP----------------------------------------------------------------
#Netdis using the Geometric-Poisson approximation as a way to obtain background expectations. 
netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = NULL)

## ----netdiszero---------------------------------------------------------------
#Netdis using no expectations (or equivalently, expectation equal to zero).
netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = 0)


