% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_net_dis.R
\name{netdis_centred_graphlet_counts}
\alias{netdis_centred_graphlet_counts}
\title{netdis_centred_graphlet_counts}
\usage{
netdis_centred_graphlet_counts(
  graphlet_counts,
  ref_ego_density_bins,
  ref_binned_graphlet_counts,
  binning_fn,
  bin_counts_fn,
  exp_counts_fn,
  max_graphlet_size
)
}
\arguments{
\item{graphlet_counts}{Ego network graphlet counts for a query graph}

\item{ref_ego_density_bins}{Either a list of previously calculated ego
network density bin edges from a reference network, or \code{NULL}, in
which case density bins are generated using the query graph itself.}

\item{ref_binned_graphlet_counts}{Either expected graphlet counts for each
ego network density bin from a reference network (a matrix with columns
labelled by graphlet ID and rows by density bin index), \code{NULL}, in
which case density binned counts are generated using the query graph itself,
or a constant numeric value to subtract from all graphlet counts.}

\item{binning_fn}{Function used to bin ego network densities. Only needed if
\code{ref_ego_density_bins} and \code{ref_binned_graphlet_counts} are
\code{NULL}. Takes densities as its single argument, and returns a named list
including keys \code{breaks} (vector of bin edges) and \code{interval_indexes}
(density bin index for each ego network).}

\item{bin_counts_fn}{Function used to calculate expected graphlet counts in
each density bin. Only needed if  \code{ref_ego_density_bins} and
\code{ref_binned_graphlet_counts} are \code{NULL}. Takes
\code{graphlet_counts}, \code{interval_indexes} (bin indexes) and
\code{max_graphlet_size} as arguments.}

\item{exp_counts_fn}{Function used to map from binned reference counts to
expected counts for each graphlet in each ego network of the query graphs.
Takes \code{ego_networks}, \code{density_bin_breaks},
\code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.}

\item{max_graphlet_size}{max graphlet size to calculate centred counts for. Currently only size 4 and 5 are supported.}
}
\value{
graphlet_counts minus exp_graphlet_counts for graphlets up to size
max_graphlet_size.
}
\description{
Calculate expected graphlet counts for each ego network in a query graph and
centre the actual counts by subtracting those calculated expected count
values.
}
