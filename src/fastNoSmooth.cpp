// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
using namespace Rcpp;

//' @title
//' Compute EMD
////'
////' @param loc1 numeric vector.
////' @param val1 numeric vector.
////' @param loc2 numeric vector.
////' @param val2 numeric vector.
//'
//' @export
// [[Rcpp::export]]
double NetEmdConstant(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
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
            temp1=(loc1[i]-curPos)*std::abs(curVal1-curVal2);
            res+=temp1;
            curVal1=val1[i];
            curPos=loc1[i];
            i+=1;
        }
        else
        {
            temp1=(loc2[j]-curPos)*std::abs(curVal1-curVal2);
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

// this can be speed up massively
std::vector<double> makeAllOffsets(NumericVector  loc1,NumericVector   loc2)
{
    int i,j;
    std::unordered_set<double> offsets1;
    for (i=0;i<loc2.size();i++)
    {
        for (j=0;j<loc1.size();j++)
        {
            offsets1.insert(loc2[i]-loc1[j]);
        }
    }
    // probably is better way to do this.
    std::vector<double> offsets(offsets1.begin(),offsets1.end());
    std::sort(offsets.begin(),offsets.end());
    return offsets;
}

//' @title
//' Compute full EMD
////'
////' @param loc1 numeric vector.
////' @param val1 numeric vector.
////' @param loc2 numeric vector.
////' @param val2 numeric vector.
//'
//' @export
// [[Rcpp::export]]
NumericVector  NetEmdFullExhaustiveNoSmooth(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
    int i,j;
    // make all of the offsets
    std::vector<double> offsets;
    offsets=makeAllOffsets(loc1,loc2);
    // Math proof that NetEMD must be less than pi
    double curBest;
    double curOffset=0;
    double curEMD;
    NumericVector loc0;
    loc0=loc1+offsets[0];
    curBest=NetEmdConstant(loc0,val1,loc2,val2);
    double bestOffset=offsets[0];
    for (i=1;i<offsets.size();i++)
    {
        loc0=loc1+offsets[i];
        curEMD=NetEmdConstant(loc0,val1,loc2,val2);
        if (curEMD<curBest)
        {
            curBest=curEMD;
            bestOffset=offsets[i];
        }
    }
    NumericVector result(2);
    result(0)=curBest;
    result(1)=bestOffset;
    return result;
}
