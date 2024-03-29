% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_net_dis.R
\name{density_binned_counts_gp}
\alias{density_binned_counts_gp}
\title{Calculate expected counts in density bins using the
geometric poisson (Polya-Aeppli) approximation.}
\usage{
density_binned_counts_gp(
  graphlet_counts,
  density_interval_indexes,
  max_graphlet_size
)
}
\arguments{
\item{graphlet_counts}{Graphlet counts for a number of ego_networks.}

\item{density_interval_indexes}{Density bin index for
each ego network.}

\item{max_graphlet_size}{Determines the maximum size of graphlets. Currently only size 4 and 5 are supported.
included in graphlet_counts.}
}
\description{
Calculate expected counts in density bins using the
geometric poisson (Polya-Aeppli) approximation.
}
