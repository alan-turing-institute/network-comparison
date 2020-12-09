
#' Netdis comparisons between one graph and many other graphs.
#'
#' Provides a heatmap and dendrogram for the network comparisons via \code{pheatmap}.
#' 
#' @param netdislist Default output of \code{netdis_many_to_many}.
#'
#' @param whatrow Selection of the row in \code{netdis_many_to_many$comp_spec} to be used for plotting.
#'
#' @param clustering_method Clusterng mehod as allowed in the \code{pheatmap} function from the \code{pheatmap} package.
#'
#' @param main Title of the plot.
#' 
#' @return Provides a heatmap and dendrogram for the network comparisons via \code{pheatmap}.
#' @export

netdis.plot <- function(netdislist,whatrow=c(1,2)[2],clustering_method="ward.D",main="Nedis"){
  attrname <- rownames(netdislist$netdis)[whatrow]
  g_NETDIScomp <- igraph::graph.edgelist(el = as.matrix(netdislist$comp_spec[,1:2]),directed = FALSE)
  edge_attr(graph = g_NETDIScomp,name = attrname) <- netdislist$netdis[whatrow,]
  adjmat <- get.adjacency(graph = g_NETDIScomp,type = "both",attr = attrname,names = TRUE,sparse = FALSE)
  vnames <- rownames(adjmat)
  
  legend1 <- seq(min(adjmat),max(adjmat),length.out = 5)
  levels1 <- round(legend1,digits = 2)
  pheatmap::pheatmap(mat = as.dist(adjmat),cluster_rows = TRUE,clustering_method = clustering_method,angle_col=45,main = main,treeheight_row = 80,labels_row = vnames,labels_col = vnames,display_numbers = TRUE,legend_breaks = legend1,legend_labels = levels1)
}
