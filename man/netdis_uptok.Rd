% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_net_dis.R
\name{netdis_uptok}
\alias{netdis_uptok}
\title{Netdis - for all graphlet sizes up to max_graphlet_size}
\usage{
netdis_uptok(
  centred_graphlet_count_vector_1,
  centred_graphlet_count_vector_2,
  max_graphlet_size
)
}
\arguments{
\item{centred_graphlet_count_vector_1}{Centred Graphlet Counts  vector for graph 1}

\item{centred_graphlet_count_vector_2}{Centred Graphlet Counts vector for graph 2}

\item{max_graphlet_size}{max graphlet size to calculate Netdis for.
The size of a graphlet is the number of nodes it contains. Netdis is
calculated for all graphlets from size 3 to size max_graphlet_size. Currently only 4 and 5 are supported.}
}
\value{
Netdis statistic calculated using centred counts for graphlets of
the specified size
}
\description{
Calculate Netdis statistic between two graphs from their Centred Graphlet
Counts (generated using \code{netdis_centred_graphlet_counts}) for all
graphlet sizes up to \code{max_graphlet_size}.
}
