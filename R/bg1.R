

getDhist <- function(x)
{
  x1=x[1:(length(x)/2)]
  x2=x[(length(x)/2+1):length(x)]
#  x1 = cumsum(x1)
#  x2 = cumsum(x2)
  
  dhist1 = dhist(x1,rep(1,length(x1)))
  dhist2 = dhist(x2,rep(1,length(x2)))
  return(list(dhist1,dhist2))
}

temp0 <- function(x)
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
  
  return(list(list(t1,t0,t2,t3),list(t1_norm,t0_norm,t2_norm,t3_norm)))
}


temp1 <- function(x)
{
  x1=x[1:(length(x)/2)]
  x2=x[(length(x)/2+1):length(x)]
  x1 = cumsum(x1)
  x2 = cumsum(x2)
  dhist1 = dhist(x1,rep(1,length(x1)))
  dhist2 = dhist(x2,rep(1,length(x2)))
  
  dhist1 = normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2 = normalise_dhist_variance(normalise_dhist_mass(dhist2))
  -min_emd(dhist1,dhist2,'optimiseRonly')$min_emd
}

temp2 <- function(n)
{
  t1=0
  print(n)
  for (i in 1:100)
  {
    start = rnorm(2*n)
    t2 = optim(start,temp1)$value
    if (t2<t1)
    {
       t1=t2
       print(t1)
    }
    
  }
  t1
}