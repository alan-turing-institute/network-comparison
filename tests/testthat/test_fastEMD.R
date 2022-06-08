context("Name of test context")

library(Rcpp)

makeFunction <- function(x1,v1,w1,x2,v2,w2)
{
  f1 <- function(x)
  {
    res1 = rep(0,length(x))
    m1 = x1[1]<x
    m2 = x<=(x1[1]+w1)
    m3 = m1*m2
    res1 = res1 + m3*(x-x1[1])*(v1[1]-0)/w1
    if (length(x1)>1)
    {
      m1 = (x1[1]+w1)<x
      m2 = x<=(x1[2])
      m3 = m1*m2
      res1 = res1 + m3*v1[1] 
      for (i in 2:length(x1))
      {
        m1 = x1[i]<x
        m2 = x<=(x1[i]+w1)
        m3 = m1*m2
        res1 = res1 + m3*((x-x1[i])*(v1[i]-v1[i-1])/w1+v1[i-1])
        if (i<length(x1)) 
        {
          m1 = (x1[i]+w1)<x
          m2 = x<=(x1[i+1])
          m3 = m1*m2
          res1 = res1 + m3*v1[i]
        }
  
      }
      
    }
    m1 = x>(x1[length(v1)]+w1)
    res1 = res1 + m1*v1[length(v1)]

    
    
    res2= rep(0,length(x))
    m1 = x2[1]<x
    m2 = x<=(x2[1]+w2)
    m3 = m1*m2
    res2 = m3*(x-x2[1])*(v2[1]-0)/w2
    
    if (length(x2)>1)
    {
      m1 = (x2[1]+w2)<x
      m2 = x<=(x2[2])
      m3 = m1*m2
      res2 = res2 + m3*v2[1] 
    }
    
    for (i in 2:length(x2))
    {
      m1 = x2[i]<x
      m2 = x<=(x2[i]+w2)
      m3 = m1*m2
      res2 = res2 + m3*((x-x2[i])*(v2[i]-v2[i-1])/w2+v2[i-1])
      if (i<length(x2)) 
      {
        m1 = (x2[i]+w2)<x
        m2 = x<=(x2[i+1])
        m3 = m1*m2
        res2 = res2 + m3*v2[i]
      }
    }
    m1 = x>(w2+x2[length(v2)])
    res2 = res2 + m1*v2[length(v2)]
    
    abs(res1-res2)
  }
  f1
}

getVal <- function(x1,v1,w1,x2,v2,w2)
{
  f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
  mi123 = min(min(x1),min(x2))
  ma123 = max(max(x1)+w1,max(x2)+w2)
  integrate(f1,mi123,ma123)
}


test_that("3 element test", {
  
            for (w1 in (1:10)/10.0)
            {
              for (w2 in (1:10)/10.0)
              {
                x1 <- c(1,2,3)
                v1 <- c(0.25,0.70,1.00)
                x2 <- c(1,2,3)
                v2 <- c(0.25,0.70,1.00)
                f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
                
                res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                
                expect_lt(abs(res1-res2),10**(-3))
              }
            }
            })

test_that("2 element test w1=0.1, w2=0.2", {
  
  w1 = 0.1
  w2 = 0.2
  x1 <- c(1,2)
  v1 <- c(0.25,0.75)
  x2 <- c(1,2)
  v2 <- c(0.5,1.00)
  f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
  res2 <- integrate(f1,0,max(x2[2],x1[2])+max(w1,w2),abs.tol=0.000000001)[[1]]
  res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
  
  expect_lt(abs(res1-res2),10**(-4))
})


test_that("1 element at 0 vs many test Mixture", {
  
  for (w1 in (1:10)*2)
  {
      x1 <- c(0)
      v1 <- c(1.00)
      x2 <- 1:w1
      v2 <- (1:w1)/w1
      f1 <- makeFunction(x1,v1,1,x2,v2,1)
      res2 <- integrate(f1,0,w1+1,abs.tol=0.000000001)[[1]]
      
      res1 <- NetEmdSmoothV2(x1,v1,1,x2,v2,1)
      
      expect_lt(abs(res1-res2),10**(-3))
    }
})


test_that("1 element vs many test Mixture", {
  
  for (w1 in (1:10)*2)
  {
      x1 <- c(w1/2)
      v1 <- c(1.00)
      x2 <- 1:w1
      v2 <- (1:w1)/w1
      f1 <- makeFunction(x1,v1,1,x2,v2,1)
      res2 <- integrate(f1,0,w1+1,abs.tol=0.000000001)[[1]]
      
      res1 <- NetEmdSmoothV2(x1,v1,1,x2,v2,1)
      
      expect_lt(abs(res1-res2),10**(-3))
    }
})



test_that("3 element test Mixture", { 
            for (w1 in (1:10)/10.0)
            {
              for (w2 in (1:10)/10.0)
              {
                x1 <- c(1,2,3)
                v1 <- c(0.65,0.70,1.00)
                x2 <- c(1,2,3)
                v2 <- c(0.25,0.70,1.00)
                f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
                
                res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                
                expect_lt(abs(res1-res2),10**(-3))
              }
            }
            })

test_that("3 element test Mixture MidPoint", { 
            w1 = 1
            w2 = 1
            for (v1_2 in (1:10)/10.0)
            {
              for (v2_2 in (1:10)/10.0)
              {
                x1 <- c(1,2,3)
                v1 <- c(0.1,v1_2,1.00)
                x2 <- c(1,2,3)
                v2 <- c(0.1,v2_2,1.00)
                f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
                
                res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                
                expect_lt(abs(res1-res2),10**(-3))
              }
            }
            })


test_that("3 element test Mixture StartPoint", { 
            w1 = 1
            w2 = 1
            for (v1_1 in (1:5)/10.0)
            {
              for (v2_1 in (1:5)/10.0)
              {
                x1 <- c(1,2,3)
                v1 <- c(v1_1,0.5,1.00)
                x2 <- c(1,2,3)
                v2 <- c(v2_1,0.5,1.00)
                f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
                
                res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                
                expect_lt(abs(res1-res2),10**(-3))
              }
            }
            })


test_that("3 element test Mixture StartLoc", { 
            w1 = 1
            w2 = 1
            for (x1_1 in (1:9)/10.0)
            {
              for (x2_1 in (1:9)/10.0)
              {
                x1 <- c(x1_1,2,3)
                v1 <- c(0.25,0.5,1.00)
                x2 <- c(x2_1,2,4)
                v2 <- c(0.3,0.5,1.00)
                f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
                
                res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                
                expect_lt(abs(res1-res2),10**(-3))
              }
            }
            })


test_that("many element test Mixture ", { 
            w1 = 1
            w2 = 1
            for (i in (2:10)*1)
            {
              for (j in (2:10)*1)
              {
                for (x123 in (1:2))
                {
                  for (y123 in (1:2))
                  {
                    x1 <- cumsum(abs(rnorm(i)))
                    v1 <- cumsum(abs(rnorm(i)))
                    w1 <- min(diff(x1))/x123
                    v1 = v1/v1[length(v1)]
                    x2 <- cumsum(abs(rnorm(j)))
                    w2 <- min(diff(x2))/y123
                    v2 <- cumsum(abs(rnorm(j)))
                    v2 = v2/v2[length(v2)]
                    f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                    top1 = max(x2[length(x2)],x1[length(x1)])+max(w1,w2)
                    bottom1 = min(x2[1],x1[1])
                    res2 <- 0 
                    res2 <- res2 + integrate(f1,bottom1,top1,abs.tol=0.000000001,subdivisions = 100000000)[[1]]
                    
                    res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                    #  if (abs(res1-res2)>0.001)
                    #  {
                    #    browser()
                    #  }
                    # Swapped to percentage error
                    expect_lt(abs(res1-res2),10**(-3))
                    
                  }
                }
              }
            }
            })


test_that("3 element test w1=0.1, w2=0.2", {
  
  w1 = 0.1
  w2 = 0.2
              x1 <- c(1,2,3)
              v1 <- c(0.25,0.70,1.00)
              x2 <- c(1,2,3)
              v2 <- c(0.25,0.70,1.00)
              f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
              res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.0000000001)[[1]]
              res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
              
              expect_lt(abs(res1-res2),10**(-4))
            })

test_that("Old failure case", {
  d1 = list()
  attr(d1,'class') <- "dhist"
  d1$locations <- 0 
  d1$masses <- 1000 
  d2 = list()
  attr(d2,'class') <- "dhist"
  d2$locations <- c(0,1,2,3) 
  d2$masses <- c(8634,1242,114,10) 
  sq1 <- netemd_single_pair(d1,d2,method='optimise',smoothing_window_width = 1)
  sq2 <- netemd_single_pair(d1,d2,method='optimiseRonly',smoothing_window_width = 1)
  expect_lt(abs(sq1$min_emd-sq2$min_emd),10**(-4))
})


test_that("Old failure case 2", {
  x1 = c(0.2862492, 0.6917626)
  v1 = c(0.6519357, 1.0000000)
  w1 = 0.2027567
  x2 = c(0.9990626, 1.4882579)
  v2 = c(0.6519357, 1.0000000)
  w2 =  0.2445976
  f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
  res1 <- integrate(f1,min(min(x1),min(x2)),max(max(x1),max(x2))+max(w1,w2))[[1]]
  res2 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
  expect_lt(abs(res2-res1),10**(-4))
})
  


test_that("Old Failure Case 2 reverse", {
  x2 = c(0.2862492, 0.6917626)
  v2 = c(0.6519357, 1.0000000)
  w2 = 0.2027567
  x1 = c(0.9990626, 1.4882579)
  v1 = c(0.6519357, 1.0000000)
  w1 =  0.2445976
  f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
  res1 <- integrate(f1,min(min(x1),min(x2)),max(max(x1),max(x2))+max(w1,w2))[[1]]
  res2 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
  expect_lt(abs(res2-res1),10**(-4))
})