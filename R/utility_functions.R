# VECTOR FUNCTIONS
rotl_vec <- function(vec, lshift) {
  num_els <- length(vec)
  select_mask <- ((1:num_els + lshift) %% num_els)
  select_mask[select_mask == 0] <- num_els
  return(vec[select_mask])
}

rotr_vec <- function(vec, rshift) {
  return(rotl_vec(vec, -rshift))
}
