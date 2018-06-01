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
   double loc1Start;
   double loc2Start;
   double loc1End;
   double loc2End;
   double h;
   res=0;
   int count123=0;
   while (1)
   {
        // lets compute the area for these segments
        if (loc1SegValStart<loc2SegValStart)
        {
            curStartVal=loc2SegValStart;
            loc2Start=loc2SegStart;
            loc1Start=loc1SegStart+(loc1SegEnd-loc1SegStart)*(loc2SegValStart-loc1SegValStart)/(loc1SegValEnd-loc1SegValStart);
        }
        else
        {
            curStartVal=loc1SegValStart;
            loc1Start=loc1SegStart;
            loc2Start=loc2SegStart+(loc2SegEnd-loc2SegStart)*(loc1SegValStart-loc2SegValStart)/(loc2SegValEnd-loc2SegValStart);
        }
        if (loc1SegValEnd<loc2SegValEnd)
        {
            curEndVal=loc1SegValEnd;
            loc1End=loc1SegEnd;
            loc2End=loc2SegStart+(loc2SegEnd-loc2SegStart)*(loc1SegValEnd-loc2SegValStart)/(loc2SegValEnd-loc2SegValStart);
        }
        else
        {
            curEndVal=loc2SegValEnd;
            loc2End=loc2SegEnd;
            loc1End=loc1SegStart+(loc1SegEnd-loc1SegStart)*(loc2SegValEnd-loc1SegValStart)/(loc1SegValEnd-loc1SegValStart);
        }
        h=(curEndVal-curStartVal)/2.0;
        if (loc1Start<loc2Start)
        {
            //case1 they is no overlap
            if (loc1End<=loc2End)
            {
                res+=(h)*(loc2Start+loc2End-loc1Start-loc1End);
            }
            else // we have a bowtie
            {
              res+=(h)*((loc2Start-loc1Start)*(loc2Start-loc1Start)+(loc1End-loc2End)*(loc1End-loc2End))/(loc1End-loc2End+loc2Start-loc1Start);
            }
        }
        else
        {
            //case1 they is no overlap
            if (loc2End<=loc1End)
            {
                res+=(h)*(loc1Start+loc1End-loc2Start-loc2End);
            }
            else // we have a bowtie
            {
                res+=(h)*((loc1Start-loc2Start)*(loc1Start-loc2Start)+(loc2End-loc1End)*(loc2End-loc1End))/(loc2End-loc1End+loc1Start-loc2Start);
            }
        }
        // update the segment under consideration
        if (i==val1.size()-1)
        {
            if (j==val2.size()-1)
            {
            //    std::cout << " i hit the break\n";
                break;
            }
            else
            {
               j+=1;
               loc2SegStart=loc2[j];
               loc2SegEnd=loc2[j]+binWidth2;
               loc2SegValStart=loc2SegValEnd;
               loc2SegValEnd=val2[j];
            }
        }
        else if ((j==val2.size()-1) || (val1[i+1]<val2[j+1]))
        {
          i+=1;
           loc1SegStart=loc1[i];
           loc1SegEnd=loc1[i]+binWidth1;
           loc1SegValStart=loc1SegValEnd;
           loc1SegValEnd=val1[i];
        }
        else
        {
           j+=1;
           loc2SegStart=loc2[j];
           loc2SegEnd=loc2[j]+binWidth2;
           loc2SegValStart=loc2SegValEnd;
           loc2SegValEnd=val2[j];
        }
   }
    return res;
}

