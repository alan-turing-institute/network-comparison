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
    m1 = x>(x1[length(v1)]+w1)
    res1 = res1 + m1*v1[length(v1)]

    
    
    res2= rep(0,length(x))
    m1 = x2[1]<x
    m2 = x<=(x2[1]+w2)
    m3 = m1*m2
    res2 = m3*(x-x2[1])*(v2[1]-0)/w2
    
    
    m1 = (x2[1]+w2)<x
    m2 = x<=(x2[2])
    m3 = m1*m2
    res2 = res2 + m3*v2[1] 
    
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

test_that("3 element test", {
  
            for (w1 in (1:10)/10.0)
            {
              for (w2 in (1:10)/10.0)
              {
                x1 <- c(1,2,3)
                v1 <- c(0.25,0.70,1.00)
                x2 <- c(1,2,3)
                v2 <- c(0.25,0.70,1.00)
                print(x1)
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
  print(x1)
  f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
  res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
  res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
  
  expect_lt(abs(res1-res2),10**(-4))
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
                print(x1)
                f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
                res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
                
                res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
                
                expect_lt(abs(res1-res2),10**(-3))
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
              print(x1)
              f1 <- makeFunction(x1,v1,w1,x2,v2,w2)
              res2 <- integrate(f1,0,max(x2[3],x1[3])+max(w1,w2),abs.tol=0.000000001)[[1]]
              res1 <- NetEmdSmoothV2(x1,v1,w1,x2,v2,w2)
              
              expect_lt(abs(res1-res2),10**(-4))
            })
