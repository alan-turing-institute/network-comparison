% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_net_dis.R
\name{ego_network_density}
\alias{ego_network_density}
\title{Calculate ego network edge densities.}
\usage{
ego_network_density(graphlet_counts)
}
\arguments{
\item{graphlet_counts}{Matrix of pre-generated graphlet, edge and node counts
(columns) for each ego network (rows). Columns must include "N" (node counts)
and "G0" (edge counts).}
}
\description{
Calculate ego network edge densities.
}
