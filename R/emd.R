#' Minimum Earth Mover's Distance (EMD)
#'
#' Calculates the minimum Earth Mover's Distance (EMD) between two discrete
#' histograms. This is the minimum EMD between the two histograms across all
#' possible offsets of histogram 1 against histogram 2.
#' @param dhist1 A \code{dhist} discrete histogram object
#' @param dhist2 A \code{dhist} discrete histogram object
#' @param method The method to use to find the minimum EMD across all potential
#' offsets for each pair of histograms. Default is "optimise" to use
#' R's built-in \code{stats::optimise} method to efficiently find the offset
#' with the minimal EMD. However, this is not guaranteed to find the global
#' minimum if multiple local minima EMDs exist. You can alternatively specify
#' the "exhaustive" method, which will exhaustively evaluate the EMD between the
#' histograms at all offsets that are candidates for the minimal EMD.
#' @return Earth Mover's Distance between the two discrete histograms
#' @export
min_emd <- function(dhist1, dhist2, method = "optimise") {
  # Require input to be a pair of "dhist" discrete histograms
  if (!(is_dhist(dhist1) && is_dhist(dhist2))) {
    stop("All inputs must be 'dhist' discrete histogram objects")
  }
  if (method == "optimise") {
    return(min_emd_optimise_fast(dhist1, dhist2))
  } else if (method == "optimiseRonly") {
    return(min_emd_optimise(dhist1, dhist2))
  } else if (method == "exhaustive") {
    return(min_emd_exhaustive(dhist1, dhist2))
  } else {
    stop("Method not recognised. Must be 'exhaustive' or ' optimise'")
  }
}



#' Minimum Earth Mover's Distance (EMD) using fast optimiser search
#'
#' Calculates the minimum Earth Mover's Distance (EMD) between two discrete
#' histograms by minimising the offset parameter of the \code{emd} function
#' using the built-in \code{stats::optimise} method.
#' @param dhist1 A \code{dhist} discrete histogram object
#' @param dhist2 A \code{dhist} discrete histogram object
#' @return Earth Mover's Distance between the two discrete histograms
#' @export
min_emd_optimise_fast <- function(dhist1, dhist2) {
  # Can we run the C++ fast implementation (only works with no smoothing)?
  if ((dhist1$smoothing_window_width == 0) &&
    (dhist2$smoothing_window_width == 0)) {
    # Determine minimum and maximum offset of range in which histograms overlap
    # (based on sliding histogram 1)
    min_offset <- min(dhist2$locations) - max(dhist1$locations)
    max_offset <- max(dhist2$locations) - min(dhist1$locations)
    # Set lower and upper range for optimise algorithm to be somewhat wider than
    # range defined by the minimum and maximum offset. This guards against a
    # couple of issues that arise if the optimise range is exactly min_offset
    # to max_offset
    # 1) If lower and upper are equal, the optimise method will throw an error
    # 2) It seems that optimise is not guaranteed to explore its lower and upper
    #    bounds, even in the case where one of them is the offset with minimum
    #    EMD
    buffer <- 0.1
    min_offset <- min_offset - buffer
    max_offset <- max_offset + buffer
    # Define a single parameter function to minimise emd as a function of offset
    val1 <- cumsum(dhist1$masses)
    val2 <- cumsum(dhist2$masses)
    val1 <- val1 / val1[length(val1)]
    val2 <- val2 / val2[length(val2)]
    loc1 <- dhist1$locations
    loc2 <- dhist2$locations
    emd_offset <- function(offset) {
      temp1 <- emd_fast_no_smoothing(loc1 + offset, val1, loc2, val2)
      temp1
    }
    # Get solution from optimiser
    soln <- stats::optimise(emd_offset,
      lower = min_offset, upper = max_offset,
      tol = .Machine$double.eps * 1000
    )
    # Return mnimum EMD and associated offset
    min_emd <- soln$objective
    min_offset <- soln$minimum
    return(list(min_emd = min_emd, min_offset = min_offset))
  } else
  {
    val1 <- cumsum(dhist1$masses)
    val2 <- cumsum(dhist2$masses)
    val1 <- val1 / val1[length(val1)]
    val2 <- val2 / val2[length(val2)]
    loc1 <- dhist1$locations
    loc2 <- dhist2$locations
    bin_width_1 <- dhist1$smoothing_window_width
    bin_width_2 <- dhist2$smoothing_window_width
    # Offset the histograms to make the alignments work
    loc1_mod <- loc1 - bin_width_1 / 2
    loc2_mod <- loc2 - bin_width_2 / 2
    # Determine minimum and maximum offset of range in which histograms overlap
    # (based on sliding histogram 1)
    min_offset <- min(loc2_mod) - max(loc1_mod) - max(bin_width_1, bin_width_2)
    max_offset <- max(loc2_mod) - min(loc1_mod) + max(bin_width_1, bin_width_2)
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
    # Define a single parameter function to minimise emd as a function of offset
    emd_offset <- function(offset) {
      temp1 <- netemd_smooth(
        loc1_mod + offset,
        val1,
        bin_width_1,
        loc2_mod,
        val2,
        bin_width_2
      )
      temp1
    }
    # Get solution from optimiser
    soln <- stats::optimise(emd_offset,
      lower = min_offset, upper = max_offset,
      tol = .Machine$double.eps * 1000
    )
    # Return mnimum EMD and associated offset
    min_emd <- soln$objective
    min_offset <- soln$minimum
    return(list(min_emd = min_emd, min_offset = min_offset))
  }
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
min_emd_optimise <- function(dhist1, dhist2) {
  # Determine minimum and maximum offset of range in which histograms overlap
  # (based on sliding histogram 1)
  min_offset <- min(dhist2$locations) - max(dhist1$locations)
  max_offset <- max(dhist2$locations) - min(dhist1$locations)

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

  # Define a single parameter function to minimise emd as a function of offset
  emd_offset <- function(offset) {
    # Construct ECMFs for each normalised histogram
    ecmf1 <- dhist_ecmf(shift_dhist(dhist1, offset))
    ecmf2 <- dhist_ecmf(dhist2)
    area_between_dhist_ecmfs(ecmf1, ecmf2)
  }

  # Get solution from optimiser
  soln <- stats::optimise(emd_offset,
    lower = min_offset, upper = max_offset,
    tol = .Machine$double.eps * 1000
  )

  # Return mnimum EMD and associated offset
  min_emd <- soln$objective
  min_offset <- soln$minimum
  return(list(min_emd = min_emd, min_offset = min_offset))
}

#' Minimum Earth Mover's Distance (EMD) using exhaustive search
#'
#' Calculates the minimum Earth Mover's Distance (EMD) between two discrete
#' histograms using an exhaustive search.
#'
#' When "sliding" two piecewise-linear empirical cumulative mass functions
#' (ECMFs) across each other to minimise the EMD between them, it is sufficient
#' to calculate the EMD at all offsets where any knots from the two ECMFs align
#' to ensure that the offset with the global minimum EMD is found.
#'
#' This is because of the piece-wise linear nature of the two ECMFs. Between any
#' two offsets where knots from the two ECMFs align, EMD will be either
#' constant, or uniformly increasing or decreasing. Therefore, there the EMD
#' between two sets of aligned knots cannot be smaller than the EMD at one or
#' other of the bounding offsets.
#' @param dhist1 A \code{dhist} discrete histogram object
#' @param dhist2 A \code{dhist} discrete histogram object
#' @return Earth Mover's Distance between the two discrete histograms
#' @export
min_emd_exhaustive <- function(dhist1, dhist2) {
  # Determine initial offset for histogram1, based on sliding histogram1
  # over histogram 2 from left to right
  step_shift <- min(dhist2$locations) - max(dhist1$locations)
  # Set initial values of minimum EMD and associated offset to impossible values
  min_emd <- Inf # Inf so that first EMD is less than this
  min_offset <- Inf
  cur_offset <- 0 # 0 so that adding first step shift gives initial offset
  # Set state variables
  distance_matrix <- NULL
  while (step_shift < Inf) {
    dhist1 <- shift_dhist(dhist1, step_shift)
    cur_offset <- cur_offset + step_shift
    cur_emd <- emd(dhist1, dhist2)
    if (cur_emd < min_emd) {
      min_emd <- cur_emd
      min_offset <- cur_offset
    }
    res <- shift_to_next_alignment(dhist1$locations, dhist2$locations,
      distance_matrix_prev = distance_matrix,
      shift_prev = step_shift
    )
    step_shift <- res$shift
    distance_matrix <- res$distance_matrix
  }
  return(list(min_emd = min_emd, min_offset = min_offset))
}

#' Minimum shift to next alignment of two location vectors
#'
#' Calculate minimum right shift of first location vector to make any pair of
#' locations from the two vectors equal
#' @param x1 First location vector. This vector is being shifted rightwards
#' @param x2 Second location vector. This vector is remaining unchanged.
#' @return Minimum non-zero right-shift to apply to x1 to align at least one
#' element of x1 with at least one element of x2
shift_to_next_alignment <- function(x1, x2, distance_matrix_prev = NULL,
                                    shift_prev = NULL) {
  if (!is.null(distance_matrix_prev) && !is.null(shift_prev)) {
    # If both distance matrix and shift from previous step provided, use these
    # to more efficiently calculate distance matrix
    distance_matrix <- (distance_matrix_prev - shift_prev)
  } else {
    # Otherwise calculate distance matrix from scratch by calculating the
    # distance from each x1 to each x2
    # NOTE: outer() generates a matrix with the first vector mapped to rows and
    # the second vector mapped to columns, so the rows will be x2 and the
    # columns x1
    distance_matrix <- outer(x2, x1, "-")
  }
  # Calculate the distance from each x1 to each x2
  # outer() generates a matrix with the first vector mapped to rows and the
  # second vector mapped to columns
  # We're stepping x1 from left to right across x2, so drop all negative
  # distances. Also drop zero distances as we want to step to the next alingment
  # even when x1 and x2 are already aligned
  distance_matrix[distance_matrix <= 0] <- Inf
  # Return the minimum positive distance as the smallest step required to align
  # any x1 with any x2
  return(list(shift = min(distance_matrix), distance_matrix = distance_matrix))
}

#' Earth Mover's Distance (EMD)
#'
#' Calculates the Earth Mover's Distance (EMD) between two discrete histograms
#' @param dhist1 A \code{dhist} discrete histogram object
#' @param dhist2 A \code{dhist} discrete histogram object
#' @return Earth Mover's Distance between the two discrete histograms
#' @export
emd <- function(dhist1, dhist2) {
  # Require the inputs to be "dhist" objects
  if (!(is_dhist(dhist1) && is_dhist(dhist2))) {
    stop("All inputs must be 'dhist' discrete histogram objects")
  }
  # Use efficient difference of cumulative histogram method that can also
  # handle non-integer bin masses and location differences
  emd_cs(dhist1, dhist2)
}

#' Earth Mover's Distance (EMD) using the difference of cumulative sums method
#'
#' Takes two discrete histograms and calculates the Wasserstein / Earth Mover's
#' Distance between the two histograms by summing the absolute difference
#' between the two cumulative histograms.
#' @references
#' Calculation of the Wasserstein Distance Between Probability Distributions on
#' the Line S. S. Vallender, Theory of Probability & Its Applications 1974 18:4,
#' 784-786 \url{http://dx.doi.org/10.1137/1118101}
#' @param dhist1 A discrete histogram as a \code{dhist} object
#' @param dhist2 A discrete histogram as a \code{dhist} object
#' @return Earth Mover's Distance between the two input histograms
#' @export
emd_cs <- function(dhist1, dhist2) {
  # Generate Empirical Cumulative Mass Functions (ECMFs) for each discrete
  # histogram
  ecmf1 <- dhist_ecmf(dhist1)
  ecmf2 <- dhist_ecmf(dhist2)
  # Calculate the area between the two ECMFs
  area <- area_between_dhist_ecmfs(ecmf1, ecmf2)
  return(area)
}

#' Earth Mover's Distance (EMD) using linear programming (LP)
#'
#' Takes two sets of histogram bin masses and bin centres and calculates the
#' Earth Mover's Distance between the two histograms by solving the Transport
#' Problem using linear programming.
#'
#' WARNING: Linear Programming solution will only give a correct answer if all
#' masses and distances between bin centres are integers.
#' @param bin_masses1 Bin masses for histogram 1
#' @param bin_masses2 Bin masses for histogram 2
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return Earth Mover's Distance between the two input histograms
#' @export
emd_lp <- function(bin_masses1, bin_masses2, bin_centres1, bin_centres2) {
  num_bins1 <- length(bin_masses1)
  num_bins2 <- length(bin_masses2)

  # Check inputs: All bins in each histogram must have a mass and centre, so
  # the bin_mass and bin_centre vectors for each histogram must have the same
  # length.
  if (length(bin_centres1) != num_bins1) {
    stop(
      "Number of bin masses and bin centres provided for histogram 1 must ",
      "be equal"
    )
  }
  if (length(bin_centres2) != num_bins2) {
    stop(
      "Number of bin masses and bin centres provided for histogram 2 must ",
      "be equal"
    )
  }

  # Generate cost matrix
  cost_mat <- cost_matrix(bin_centres1, bin_centres2)

  # Linear Programming solver requires all bin masses and transportation costs
  # to be integers to generate correct answer
  if (!isTRUE(all.equal(bin_masses1, floor(bin_masses1)))) {
    stop(
      "All bin masses for histogram 1 must be integers for accurate Linear ",
      "Programming solution"
    )
  }
  if (!isTRUE(all.equal(bin_masses2, floor(bin_masses2)))) {
    stop(
      "All bin masses for histogram 2 must be integers for accurate ",
      "Linear Programming solution"
    )
  }
  if (!isTRUE(all.equal(cost_mat, floor(cost_mat)))) {
    stop("All costs must be integers for accurate Linear Programming solution")
  }
  row_signs <- rep("==", num_bins1)
  col_signs <- rep("<=", num_bins2)
  s <- lpSolve::lp.transport(
    cost.mat = cost_mat, row.signs = row_signs,
    col.signs = col_signs, row.rhs = bin_masses1,
    col.rhs = bin_masses2
  )
  return(s$objval)
}

#' Inter-bin cost matrix from bin centres
#'
#' Generates a matrix for the cost of moving a unit of mass between each bin in
#' histogram 1 and each bin in histogram 2.
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return Cost matrix
cost_matrix <- function(bin_centres1, bin_centres2) {
  # Calculate distances between all bins in network 1 and all bins in network 2
  num_bins1 <- length(bin_centres1)
  num_bins2 <- length(bin_centres2)
  loc_mat1 <- matrix(
    bin_centres1,
    nrow = num_bins1,
    ncol = num_bins2,
    byrow = FALSE
  )
  loc_mat2 <- matrix(
    bin_centres2,
    nrow = num_bins1,
    ncol = num_bins2,
    byrow = TRUE
  )
  cost_mat <- abs(loc_mat1 - loc_mat2)
  return(cost_mat)
}
