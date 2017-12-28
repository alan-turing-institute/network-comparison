library('Rcpp')
sourceCpp('cppdhist.cpp')
#sourceCpp('cppConstantUpdateVersion.cpp')
#sourceCpp('constantExhaustiveCombined.cpp')

giveMeHist <-function(i)
{
    d1=graph_features_to_histograms(matrix(1:10000000,nrow=i,ncol=2))
    d1=d1[[1]]
    d1=normalise_dhist_mass(d1)
    d1=normalise_dhist_variance(d1)
    d1=shift_dhist(d1,-dhist_mean_location(d1))
d1
}

runTest <- function(lim)
{
    res=c()
    for (i in 5:lim)
    {
        print('values of i')
        print(c(i))
        d1=graph_features_to_histograms(matrix(1:10000,nrow=i,ncol=2))
        d2=graph_features_to_histograms(matrix(1:10000,nrow=i-1,ncol=2))
res1=net_emd_fast(d1,d2,method='exhaust')
res=append(res,res1)
    }
    res
}

runTest2 <- function(lim)
{
    res=c()
    for (i in 5:lim)
    {
        print('values of i')
        print(c(i))
        d1=graph_features_to_histograms(matrix(sample.int(10000,size=i*2+1),nrow=i,ncol=2))
        d2=graph_features_to_histograms(matrix(sample.int(10000,size=2*i+2),nrow=i-1,ncol=2))
res1=net_emd_fast(d1,d2,method='exhaust')
res=append(res,res1)
    }
    res
}


emdOpt <- function( x )
{
    x1=x[1:(length(x)/2)]
    x2=x[(1+length(x)/2):length(x)]
    m1=matrix(x1,nrow=length(x1),ncol=1)
    m1=cbind(m1,m1)
    m2=matrix(x2,nrow=length(x2),ncol=1)
    m2=cbind(m2,m2)


    d1=graph_features_to_histograms(m1)
    d2=graph_features_to_histograms(m2)
    d1=d1[[1]]
    d2=d2[[1]]
    d1=normalise_dhist_mass(d1)
    d1=normalise_dhist_variance(d1)
    d1=shift_dhist(d1,-dhist_mean_location(d1))
    d2=normalise_dhist_mass(d2)
    d2=normalise_dhist_variance(d2)
    d2=shift_dhist(d2,-dhist_mean_location(d2))
    res=emd(d1,d2)
    -res
}

runTest4 <- function()
{
    result=c()
    for (i in 5:100)
    {
        curBest=0
        for (rep in 1:10)
        {
            qw1=sample.int(1000,2*i)
    res=optim(qw1,emdOpt)
            if (res$value<curBest)
            {
                curBest=res$value
            }
        }
        print(c(i,-curBest))
    result=append(result,-curBest)
    }
    result
}


runTest3 <- function(lim)
{
    i=10
    result=c()
    print('values of i')
    print(c(i))
    curBest=0;
    for (i in 5:1000)
    {
        curBest=0
for (rep in 1:100)
{
    print(c(i,rep,curBest))
    d1=graph_features_to_histograms(matrix(sample.int(10000,size=i*2+1),nrow=i,ncol=2))
    d2=graph_features_to_histograms(matrix(sample.int(10000,size=2*i+2),nrow=i-1,ncol=2))
    d1=d1[[1]]
    d2=d2[[1]]
    d1=normalise_dhist_mass(d1)
    d1=normalise_dhist_variance(d1)
    d1=shift_dhist(d1,-dhist_mean_location(d1))
    d2=normalise_dhist_mass(d2)
    d2=normalise_dhist_variance(d2)
    d2=shift_dhist(d2,-dhist_mean_location(d2))
    res=emd(d1,d2)
    if (res>curBest)
    {
        curBest=res
    }
  }  
   result=append(result,curBest) 
    }
    result
}

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
  } else if(method == "exhaustiveHalf"){
    return(min_emd_exhaustive_fastHalf(dhist1, dhist2))
  } else if(method == "exhaustiveHalfNoOffset"){
    return(min_emd_exhaustive_fastHalfNoOffset(dhist1, dhist2))
  } else if(method == "exhaustiveVer2"){
    return(min_emd_exhaustive_fastVer2(dhist1, dhist2))
  } else if(method == "exhaustiveVer3"){
    return(min_emd_exhaustive_fastVer3(dhist1, dhist2))
  } else if(method == "exhaust"){
    return(min_emd_exhaustive_fastVerTest(dhist1, dhist2))
  } else {
    stop("Method not recognised. Must be 'exhaustive' or ' optimise'")
  }
}



min_emd_exhaustive_fastVerTest <- function(dhist1, dhist2) {
    print("hi")
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dhist1$locations
      loc2=dhist2$locations
      offsets=matrix(nrow=length(loc1)*length(loc2),ncol=1)
      count=1
      for (i in 1:length(loc1))
      {
          for (j in 1:length(loc2))
          {
              offsets[count]=loc2[j]-loc1[i]
              count=count+1
          }
      }
      print("bye")
      offsets=unique(offsets)
      offsets=sort(offsets)
      curBest=0.0001+emd( shift_dhist(dhist1,min(abs(offsets))),dhist2) 
      print(curBest)
      curLimit=-10000
      count=0
      l0=list()
      l1=matrix(nrow=10000,ncol=2)
      l1Count=0
      offsets=c(max(offsets),offsets)
      countSkip=0
      count0=0
      for (i in 1:length(offsets))
      {
          if (offsets[i]<curLimit)
          {
              count0=count0+1
              next
          }
          h=1
          if (l1Count>0)
          {
              for (j in 1:l1Count)
              {
                  if (offsets[i]<l1[j,2])
                  {
                      if (offsets[i]>l1[j,1])
                      {
                          h=0
                      }
                  }
              }
          }
          if (h==0)
          {
              countSkip=countSkip+1
          }
          if (h==1) #(offsets[i]>=curLimit)
          {
              countSkip=0
              count=count+1
          temp1=emd( shift_dhist(dhist1,offsets[i]),dhist2) 
              if (temp1<curBest)
              {
                  curBest=temp1
              ##    for (j in 1:length(l0))
              ##    {
              ##        t1=l0[j,1]/fg1
              ##        t2=l0[j,2]/fg1-curBest
              ##        l1[[length(l1)+1]]=c(t1-t2,t1+t2)
              ##    }
              }
              t5=temp1/curBest
              t6=offsets[i]/t5 
                  for (j in 1:10)
                  {
                      fg1=2*j+1
     #         if (temp1/fg1>curBest)
     #         {
     #                 t1=offsets[i]/fg1
     #                 t2=(temp1/fg1-curBest)
     #                 l1[l1Count+1,]=c(t1-t2,t1+t2)
     #                 l1Count=l1Count+1
     #                 if (l1Count>nrow(l1)-1)
     #                 {
     #                     l1=rbind(l1,matrix(nrow=10000,ncol=2))
     #                 }
     #         ##        l1= l1[!(l1[,1]>t1-t2 & l1[,2]<t1+t2),]
     #             }
              }
              if (i>1)
              {
              t1=offsets[i]
              t2=temp1-curBest
#              l1[l1Count+1,]=c(t1-t2,t1+t2)
              curLimit=t1+t2
    #          if (t6>curLimit)
    #          {curLimit=t6}
              }
          }
      }
      print(c(count,count0,length(offsets)))
      curBest <- (count/(count+count0))
  return(list(min_emd = count, min_offset = 1/0))
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

min_emd_exhaustive_fastHalf <- function(dhist1, dhist2) {
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dhist1$locations
      loc2=dhist2$locations
      temp1=constantVersionExhaustiveHalf(loc1,val1,loc2,val2)
  return(list(min_emd = temp1, min_offset = 1/0))
}


min_emd_exhaustive_fastHalfNoOffset <- function(dhist1, dhist2) {
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dhist1$locations
      loc2=dhist2$locations
      temp1=constantVersionExhaustiveHalfNoOffsetCalc(loc1,val1,loc2,val2)
  return(list(min_emd = temp1, min_offset = 1/0))
}


min_emd_exhaustive_fastVer2 <- function(dhist1, dhist2) {
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dput(dhist1$locations)
      loc2=dput(dhist2$locations)
      temp1=constantVersionWithUpdates(loc1,val1,loc2,val2)
  return(list(min_emd = temp1, min_offset = 1/0))
}


min_emd_exhaustive_fastVer3 <- function(dhist1, dhist2) {
      val1 <- cumsum(dhist1$masses)
      val2 <- cumsum(dhist2$masses)
      val1 <- val1/val1[length(val1)]
      val2 <- val2/val2[length(val2)]
      loc1=dput(dhist1$locations)
      loc2=dput(dhist2$locations)
      temp1=constantVersionExhaustiveCombined(loc1,val1,loc2,val2)
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
      count=0
      emd_offset <- function(offset) {
          print(offset)
          print(count)
          count=count+1
        temp1<- constantVersion(loc1+offset,val1,loc2,val2)
        temp1
      }
      soln <- stats::optimise(emd_offset, lower = min_offset, upper = max_offset, 
                              tol = .Machine$double.eps*1000)
      print("number of runs")
      print(count)
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

