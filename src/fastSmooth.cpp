// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

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
double NetEmdSmooth(NumericVector loc1,NumericVector val1,double binWidth1,NumericVector loc2,NumericVector val2,double binWidth2)
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
   // stores the result
   res=0;
   //TODO be worried about adding lots of small numbers

   // current location on hist 1 and hist 2
   i=0;
   j=0;
   double cdfLower=0;
   double loc1SegStart=loc1[0];
   double loc1SegEnd=loc1[0]+binWidth1;
   double loc1SegValStart=0;
   double loc1SegValEnd=val1[0];
   double loc2SegStart=loc2[0];
   double loc2SegEnd=loc2[0]+binWidth2;
   double loc2SegValStart=0;
   double loc2SegValEnd=val2[0];
   double curStartVal;
   double curEndVal;
   res=0;
   while (1)
   {
        // lets compute the area for these segments
        if (loc1SegValStart<loc2SegValStart)
        {
            curStartVal=loc2SegValStart;
            loc2Start=loc2SegStart
            loc1Start=loc1SegStart+(loc1SegEnd-loc1SegValStart)*(loc2SegValStart-loc1SegValStart)/(loc1SegValEnd-loc1SegValStart)
        }
        else
        {
            curStartVal=loc1SegValStart;
            loc1Start=loc1SegStart
            loc2Start=loc2SegStart+(loc2SegEnd-loc2SegValStart)*(loc1SegValStart-loc2SegValStart)/(loc2SegValEnd-loc2SegValStart)
        }
        if (loc1SegValEnd<loc2SegValEnd)
        {
            curEndVal=loc1CurValEnd;
            loc1End=loc1SegEnd
            loc2End=loc2SegStart+(loc2SegEnd-loc2SegValStart)*(loc1SegValEnd-loc2SegValStart)/(loc2SegValEnd-loc2SegValStart)
        }
        else
        {
            curEndVal=loc2CurValEnd;
            loc2End=loc2SegEnd
            loc1End=loc1SegStart+(loc1SegEnd-loc1SegValStart)*(loc2SegValEnd-loc1SegValStart)/(loc1SegValEnd-loc1SegValStart)
        }
        // okay so we now have our bounds
        h=curEndVal-curStartVal;
        if (loc1Start<loc2Start)
        {
            //case1 they is no overlap
            if (loc1End<=loc2End)
            {
                res+=(h/2.0)*(loc2Start+loc2End-loc1Start-loc1End)
            }
            else // we have a bowtie
            {
                res+=(h/2.0)*((loc2Start-loc1Start)*(loc2Start-loc1Start)+(loc1End-loc2End)*(loc1End-loc2End))/(loc1End-loc2End+loc2Start-loc1Start)
            }
        }
        else
        {
            //case1 they is no overlap
            if (loc2End<=loc1End)
            {
                res+=(h/2.0)*(loc1Start+loc1End-loc2Start-loc2End)
            }
            else // we have a bowtie
            {
                res+=(h/2.0)*((loc1Start-loc2Start)*(loc1Start-loc2Start)+(loc2End-loc1End)*(loc2End-loc1End))/(loc2End-loc1End+loc1Start-loc2Start)
            }
        }
        // update the segment under consideration
        if (i==val1.size()-1)
        {
            if (j==val2.size()-1)
            {
                break;
            }
            else
            {
               j+=1;
               loc2SegStart=loc2[j];
               loc2SegEnd=loc2[j]+binWidth2;
               loc2SegValStart=loc2CurValEnd;
               loc2SegValEnd=val2[j];
            }
        }
        else if (j==val2.size()-1)
        {
           loc1SegStart=loc1[i];
           loc1SegEnd=loc1[i]+binWidth1;
           loc1SegValStart=loc1CurValEnd;
           loc1SegValEnd=val1[i];
        }
        else if (val1[i+1]<val2[j+1])
        {
            i+=1;
            // update the i segment
           loc1SegStart=loc1[i];
           loc1SegEnd=loc1[i]+binWidth1;
           loc1SegValStart=loc1CurValEnd;
           loc1SegValEnd=val1[i];
        }
        else
        {
           j+=1;
           loc2SegStart=loc2[j];
           loc2SegEnd=loc2[j]+binWidth2;
           loc2SegValStart=loc2CurValEnd;
           loc2SegValEnd=val2[j];
        }
   }
    return res;
}

