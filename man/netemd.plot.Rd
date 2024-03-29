% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingFunctions.R
\name{netemd.plot}
\alias{netemd.plot}
\title{Heatmap of NetEmd many-to-many comparisons}
\usage{
netemd.plot(
  netemdlist,
  clustering_method = "ward.D",
  main = "NetEmd",
  docluster = TRUE
)
}
\arguments{
\item{clustering_method}{Clustering method as allowed in the \code{pheatmap} function from the \code{pheatmap} package. The dendrogram will appear if \code{docluster} is TRUE (default).}

\item{main}{Title of the plot.}

\item{docluster}{controls the order of the rows and columns. If TRUE (default) the rows and columns will be reordered to create the dendrogram. If FALSE, then only the heatmap is drawn.}

\item{netdislist}{Default output of \code{netdis_many_to_many}.}

\item{whatrow}{Selection of the row in \code{netdis_many_to_many$comp_spec} to be used for plotting.}
}
\value{
Provides a heat map and dendrogram for the network comparisons via \code{pheatmap}.
}
\description{
Provides a heatmap and dendrogram for the network comparisons via \code{pheatmap}.
}
