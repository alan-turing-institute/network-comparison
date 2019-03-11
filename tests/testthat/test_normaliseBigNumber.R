computeNetEMDs <- function(x)
{
  x1=x[1:(length(x)/2)]
  x2=x[(length(x)/2+1):length(x)]
  
  dhist1 = dhist(x1,rep(1,length(x1)))
  dhist2 = dhist(x2,rep(1,length(x2)))
  
  t3 = net_emd(dhist1,dhist2,'optimise')
  
  dhist1 = normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2 = normalise_dhist_variance(normalise_dhist_mass(dhist2))
  
  t2 = min_emd(dhist1,dhist2,'optimise')$min_emd
  
  dhist1 = dhist(x1,rep(1,length(x1)))
  dhist2 = dhist(x2,rep(1,length(x2)))
  t0 = net_emd(dhist1,dhist2,'optimiseRonly')
  
  dhist1 = normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2 = normalise_dhist_variance(normalise_dhist_mass(dhist2))
  t1 = min_emd(dhist1,dhist2,'optimiseRonly')$min_emd
  
  # Mean centred version 
  dhist1 = dhist(x1,rep(1,length(x1)))
  dhist2 = dhist(x2,rep(1,length(x2)))
  
  dhist1 = mean_centre_dhist(dhist1)
  dhist2 = mean_centre_dhist(dhist2)
  t3_norm = net_emd(dhist1,dhist2,'optimise')
  
  dhist1 = normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2 = normalise_dhist_variance(normalise_dhist_mass(dhist2))
  
  t2_norm = min_emd(dhist1,dhist2,'optimise')$min_emd
  
  dhist1 = dhist(x1,rep(1,length(x1)))
  dhist2 = dhist(x2,rep(1,length(x2)))
  
  dhist1 = mean_centre_dhist(dhist1)
  dhist2 = mean_centre_dhist(dhist2)
  t0_norm = net_emd(dhist1,dhist2,'optimiseRonly')
  
  
  dhist1 = normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2 = normalise_dhist_variance(normalise_dhist_mass(dhist2))
  t1_norm = min_emd(dhist1,dhist2,'optimiseRonly')$min_emd
  
  return(list(t1,t0,t2,t3,t1_norm,t0_norm,t2_norm,t3_norm))
}


context("Big Number Normalisation - 3 Big Numbers")
test_that("Big small and big small", {
  result <- computeNetEMDs(c(1,2,3,10000000000,-10000000000,2,3,100000000000))
  for (i in 1:8)
  {
      for (j in 1:8)
      {
        expect_true(results[[i]]==results[[j]])
      }
  }
})


context("Big Number Normalisation - 1 Big Numbers")
test_that("One big Number", {
  result <- computeNetEMDs(c(1,2,3,10000000000,1,2,3,4))
  for (i in 1:8)
  {
      for (j in 1:8)
      {
        expect_true(results[[i]]==results[[j]])
      }
  }
})

context("Big Number Normalisation - 2 Big Numbers")
test_that("Two big Numbers", {
  result <- computeNetEMDs(c(1,2,3,10000000000,100000000000,2,3,4))
  for (i in 1:8)
  {
      for (j in 1:8)
      {
        expect_true(results[[i]]==results[[j]])
      }
  }
})