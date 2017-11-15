// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
using namespace Rcpp;

double abs(double x)
{
    if (x>0)
    {return x;}
    else
    {return -x;}
}
// [[Rcpp::export]]
double constantVersionExhaustive(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
    // this is an experiment if this will be faster than running the standard exhaustive search
    std::unordered_set<double> offsets1;
    for (i=0;i<loc1.size();i++)
    {
        for (j=0;j<loc2.size();j++)
        {
            offsets1.insert(loc1[i]-loc2[j]);
        }
    }

}
// [[Rcpp::export]]
double constantVersion(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
  //init
   double res=0;
   double curVal1,curVal2;
   double curPos;
   double temp1;
   int count;
   int i,j,k;
   //place start of windows before 
   //start of histogram so we can start the loop
   if (loc1[0]<loc2[0])
   {
       curPos=loc1[0]-1.0;
   }
   else
   {
       curPos=loc2[0]-1.0;
   }
   // current value of histogram 1 an 2
   curVal1=0;
   curVal2=0;
   // stores the result
   res=0; 
   //TODO be worried about adding lots of small numbers
   
   // current location on hist 1 and hist 2
   i=0;
   j=0;
    while (1)
    {
        if (i==loc1.size())
        {break;}
        if (j==loc2.size())
        {break;}
        if (loc1[i]<loc2[j])
        {
            temp1=(loc1[i]-curPos)*abs(curVal1-curVal2);
            res+=temp1;
            curVal1=val1[i];
            curPos=loc1[i];
            i+=1;
        }
        else
        {
            temp1=(loc2[j]-curPos)*abs(curVal1-curVal2);
            res+=temp1;
            curVal2=val2[j];
            curPos=loc2[j];
            j+=1;
        }
    }
    if (i<loc1.size())
    {
        for (k=i;k<loc1.size();k++)
        {
            res+=(loc1[k]-curPos)*(1.0-curVal1);
            curVal1=val1[k];
            curPos=loc1[k];
        }
    }
    else
    {
        for (k=j;k<loc2.size();k++)
        {
            res+=(loc2[k]-curPos)*(1.0-curVal2);
            curVal2=val2[k];
            curPos=loc2[k];
        }
    }
    return res;
}

