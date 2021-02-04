
#' Heatmap of Netdis many-to-many comparisons
#'
#' Provides a heatmap and dendrogram for the network comparisons via \code{pheatmap}.
#' 
#' @param netdislist Default output of \code{netdis_many_to_many}.
#'
#' @param whatrow Selection of the row in \code{netdis_many_to_many$comp_spec} to be used for plotting.
#'
#' @param clustering_method Clustering method as allowed in the \code{pheatmap} function from the \code{pheatmap} package. The dendrogram will appear if \code{docluster} is TRUE (default).
#'
#' @param main Title of the plot.
#' 
#' @param docluster controls the order of the rows and columns. If TRUE (default) the rows and columns will be reordered to create the dendrogram. If FALSE, then only the heatmap is drawn.
#' 
#' @return Provides a heatmap and dendrogram for the network comparisons via \code{pheatmap}.
#' @export

netdis.plot <- function(netdislist,whatrow=c(1,2)[2],clustering_method="ward.D",main="Nedis",docluster=TRUE){
  adjmat <- cross_comp_to_matrix(measure = netdislist$netdis[whatrow,], cross_comparison_spec = netdislist$comp_spec)
  vnames <- rownames(adjmat)
  
  legend1 <- seq(min(adjmat),max(adjmat),length.out = 5)
  levels1 <- round(legend1,digits = 2)
  pheatmap::pheatmap(mat = as.dist(adjmat),cluster_rows = docluster,cluster_cols = docluster,clustering_method = clustering_method,angle_col=45,main = main,treeheight_row = 80,labels_row = vnames,labels_col = vnames,display_numbers = TRUE,legend_breaks = legend1,legend_labels = levels1)
}




#' Heatmap of NetEmd many-to-many comparisons
#'
#' Provides a heatmap and dendrogram for the network comparisons via \code{pheatmap}.
#' 
#' @param netdislist Default output of \code{netdis_many_to_many}.
#'
#' @param whatrow Selection of the row in \code{netdis_many_to_many$comp_spec} to be used for plotting.
#'
#' @param clustering_method Clustering method as allowed in the \code{pheatmap} function from the \code{pheatmap} package. The dendrogram will appear if \code{docluster} is TRUE (default).
#'
#' @param main Title of the plot.
#' 
#' @param docluster controls the order of the rows and columns. If TRUE (default) the rows and columns will be reordered to create the dendrogram. If FALSE, then only the heatmap is drawn.
#' 
#' @return Provides a heat map and dendrogram for the network comparisons via \code{pheatmap}.
#' @export

netemd.plot <- function(netemdlist,clustering_method="ward.D",main="NetEmd",docluster=TRUE){
  adjmat <- cross_comp_to_matrix(measure = netemdlist$netemds, cross_comparison_spec = netemdlist$comp_spec)
  vnames <- rownames(adjmat)
  
  legend1 <- seq(min(adjmat),max(adjmat),length.out = 5)
  levels1 <- round(legend1,digits = 2)
  pheatmap::pheatmap(mat = as.dist(adjmat),cluster_rows = docluster,cluster_cols = docluster,clustering_method = clustering_method,angle_col=45,main = main,treeheight_row = 80,labels_row = vnames,labels_col = vnames,display_numbers = TRUE,legend_breaks = legend1,legend_labels = levels1)
  
}
