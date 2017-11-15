library('Rcpp')
sourceCpp('cppdhist.cpp')

min_emd_fast <- function(dhist1, dhist2, method = "optimise") {
  # Require input to be a pair of "dhist" discrete histograms 
  if(!(is_dhist(dhist1) && is_dhist(dhist2))) {
    stop("All inputs must be 'dhist' discrete histogram objects")
  }

  type1 <- attr(dhist1, "type")
  type2 <- attr(dhist2, "type")
  type1 <- 'constant'
  type2 <- 'constant'
  if(type1 != type2) {
    stop("dhists must have the same type")
  }
  if(method == "optimise") {
    return(min_emd_optimise_fast(dhist1, dhist2))
  } else if(method == "exhaustive"){
    return(min_emd_exhaustive_fast(dhist1, dhist2))
  } else {
    stop("Method not recognised. Must be 'exhaustive' or ' optimise'")
  }
}



min_emd_exhaustive_fast <- function(dhist1, dhist2) {
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dhist1$locations
      loc2=dhist2$locations
      temp1=constantVersionExhaustive(loc1,val1,loc2,val2)
  return(list(min_emd = temp1, min_offset = 1/0))
}





compute_emd_offset <- function(dhist1,dhist2,offset) {
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dhist1$locations
      loc2=dhist2$locations
        temp1<- constantVersion(loc1+offset,val1,loc2,val2)
        temp1
      }


#' Minimum Earth Mover's Distance (EMD) using optimiser search
#' 
#' Calculates the minimum Earth Mover's Distance (EMD) between two discrete 
#' histograms by minimising the offset parameter of the \code{emd} function 
#' using the built-in \code{stats::optimise} method.
#' @param dhist1 A \code{dhist} discrete histogram object
#' @param dhist2 A \code{dhist} discrete histogram object
#' @return Earth Mover's Distance between the two discrete histograms
#' @export



min_emd_optimise_fast <- function(dhist1, dhist2) {
  # Determine minimum and maximum offset of range in which histograms overlap
  # (based on sliding histogram 1)
  min_offset <- min(dhist2$locations) - max(dhist1$locations)
  max_offset <- max(dhist2$locations) - min(dhist1$locations)
  type <- "constant"

  # Set lower and upper range for optimise algorithm to be somewhat wider than
  # range defined by the minimum and maximum offset. This guards against a
  # couple of issues that arise if the optimise range is exactly min_offset 
  # to max_offset
  # 1) If lower and upper are equal, the optimise method will throw and error
  # 2) It seems that optimise is not guaranteed to explore its lower and upper
  #    bounds, even in the case where one of them is the offset with minimum
  #    EMD
  buffer <- 0.1
  min_offset <- min_offset - buffer
  max_offset <- max_offset + buffer
  if(type == "constant") {
      # Define a single parameter function to minimise emd as a function of offset
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dhist1$locations
      loc2=dhist2$locations
      emd_offset <- function(offset) {
        temp1<- constantVersion(loc1+offset,val1,loc2,val2)
        temp1
      }
      soln <- stats::optimise(emd_offset, lower = min_offset, upper = max_offset, 
                              tol = .Machine$double.eps*1000)
  }
  else
  {
    stop("Type not recognised")
  }
  # Get solution from optimiser
  
  # Return mnimum EMD and associated offset
  min_emd <- soln$objective
  min_offset <- soln$minimum
  print(c(min_emd,min_offset))
  return(list(min_emd = min_emd, min_offset = min_offset))
}


net_emd_single_pair_fast <- function(dhist1, dhist2, method = "optimise", 
                                smoothing_window_width = 0) {
  # Present dhists as smoothed or unsmoothed histograms depending on the value 
  # of smoothing_window_width
  # NOTE: This MUST be done prior to any variance normalisation as the 
  # calculation of variance differs depending on whether or not the histograms 
  # are smoothed (i.e. we need to ensure that the smoothing_window_width 
  # attribute of the dhists is set to the smoothing_window_width parameter
  # provided by the caller)
  # TODO: Consider moving the smoothing of histograms outside to the user's 
  # calling code. It feels a bit untidy in here.
  if(smoothing_window_width == 0) {
    dhist1 <- as_unsmoothed_dhist(dhist1)
    dhist2 <- as_unsmoothed_dhist(dhist2)
  } else {
    dhist1 <- as_smoothed_dhist(dhist1, smoothing_window_width)
    dhist2 <- as_smoothed_dhist(dhist2, smoothing_window_width)
  }
  ## add mean centering

  ### Stores means to fix offset later
  mean1 <- dhist_mean_location(dhist1)
  mean2 <- dhist_mean_location(dhist2)
  ### Mean centre
  dhist1<-mean_centre_dhist(dhist1)
  dhist2<-mean_centre_dhist(dhist2)

  var1 <- dhist_variance(dhist1)
  var2 <- dhist_variance(dhist2)
  # Normalise histogram to unit mass and unit variance
  dhist1_norm <- normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2_norm <- normalise_dhist_variance(normalise_dhist_mass(dhist2))
  
  result <- min_emd_fast(dhist1_norm, dhist2_norm, method = method)
  result$min_offset <- result$min_offset +mean2/var2-mean1/var1
  return(result)
}



net_emd_fast <- function(dhists1, dhists2, method = "optimise", 
                    return_details = FALSE, smoothing_window_width = 0) {
  # Require either a pair of "dhist" discrete histograms or two lists of "dhist"
  # discrete histograms
  pair_of_dhist_lists <- all(purrr::map_lgl(dhists1, is_dhist)) && all(purrr::map_lgl(dhists2, is_dhist))
  
  # If input is two lists of "dhist" discrete histograms, determine the minimum
  # EMD and associated offset for pairs of histograms taken from the same 
  # position in each list
  if(pair_of_dhist_lists) {
    details <- purrr::map2(dhists1, dhists2, function(dhist1, dhist2) {
      net_emd_single_pair_fast(dhist1, dhist2, method = method,
                          smoothing_window_width = smoothing_window_width)
      })
    # Collect the minimum EMDs and associated offsets for all histogram pairs
    min_emds <- purrr::simplify(purrr::transpose(details)$min_emd)
    min_offsets <- purrr::simplify(purrr::transpose(details)$min_offset)
    # The NetEMD is the arithmetic mean of the minimum EMDs for each pair of 
    # histograms
    arithmetic_mean <- sum(min_emds) / length(min_emds)
    net_emd <- arithmetic_mean
    # Return just the NetEMD or a list including the NetEMD plus the details of
    # the minumum EMD and associated offsets for the individual histograms
    # Note that the offsets represent shifts after the histograms have been
    # scaled to unit variance
    if(return_details) {
      return(list(net_emd = net_emd, min_emds = min_emds, min_offsets = min_offsets))
    } else {
      return(arithmetic_mean)
    }
  }
  else {
    # Wrap each member of a single pair of histograms is a list and recursively
    # call this net_emd function. This ensures they are treated the same.
    return(net_emd_fast(list(dhists1), list(dhists2), method = method, 
                   return_details = return_details,
                   smoothing_window_width = smoothing_window_width))
  }
}

