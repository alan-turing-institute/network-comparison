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
   // hist1 variables
   double loc1SegStart=loc1[0];  //- start of a Segment in x
   double loc1SegEnd=loc1[0]+binWidth1; //- end of a Segment in x
   double loc1SegValStart=0; //- start of a Segment in y
   double loc1SegValEnd=val1[0]; //- end of a Segment in y

   // hist2 variables
   double loc2SegStart=loc2[0];
   double loc2SegEnd=loc2[0]+binWidth2;
   double loc2SegValStart=0;
   double loc2SegValEnd=val2[0];

   double curStartVal; // start value in y 
   double curEndVal; // end value in y 
   double loc1Start; // start value in x hist1 
   double loc2Start; // start value in x hist2
   double loc1End; // end value in x hist1 
   double loc2End; // end value in x hist2
   double h;
   res=0;
   // set as 0 as at bottom of hist
   curStartVal=0;

   // need to know if first y segment ends with hist1 or hist2 
   // Need to set the first start locations 
   // Commented this section as they are both set to zero
   if (loc1SegValStart<loc2SegValStart)
   {
     
     loc2Start=loc2SegStart;
     loc1Start=loc1SegStart+(loc1SegEnd-loc1SegStart)*(loc2SegValStart-loc1SegValStart)/(loc1SegValEnd-loc1SegValStart);
   }
   else
   {
     loc1Start=loc1SegStart;
     loc2Start=loc2SegStart+(loc2SegEnd-loc2SegStart)*(loc1SegValStart-loc2SegValStart)/(loc2SegValEnd-loc2SegValStart);
   }
   // loc1Start=loc1SegStart;
   // loc2Start=loc2SegStart;
   while (1)
   {
        // lets compute the area for this segments
        //
        // Case where hist1  ends first 
        if (loc1SegValEnd<loc2SegValEnd)
        {
            curEndVal=loc1SegValEnd; // end of current segment in y
            loc1End=loc1SegEnd; // end of current segment in x in dhist1
            // compute where end of current segment will be in x in dhist2 proportionally
            loc2End=loc2SegStart+(loc2SegEnd-loc2SegStart)*(loc1SegValEnd-loc2SegValStart)/(loc2SegValEnd-loc2SegValStart);
        }
        else // hist2 comes first
        {
            curEndVal=loc2SegValEnd; // end of current segment in y
            loc2End=loc2SegEnd; // end of current segment in x in dhist2
            // compute where end of current segment will be in x in dhist1 proportionally
            loc1End=loc1SegStart+(loc1SegEnd-loc1SegStart)*(loc2SegValEnd-loc1SegValStart)/(loc1SegValEnd-loc1SegValStart);
        }
        h=(curEndVal-curStartVal)/2.0;
        if (loc1Start<loc2Start)
        {
            //case1 they is no overlap
            if (loc1End<=loc2End)
            {
		// difference between two right triangles
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
		// difference between two right triangles
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
		// we are at the end of both segments
                break;
            }
            else
            {
	       // at the end of dhist 1 so move dhist2
               j+=1;
               loc2SegStart=loc2[j];
               loc2SegEnd=loc2[j]+binWidth2;
               loc2SegValStart=loc2SegValEnd;
               loc2SegValEnd=val2[j];
            }
        }
        else if ((j==val2.size()-1) || (val1[i+1]<val2[j+1]))
        {
	  // at the end of dhist 2 so move dhist1
	  // OR
	  // Next segment in dhist1 is smaller than the next segment in dhist2
          i+=1;
           loc1SegStart=loc1[i];
           loc1SegEnd=loc1[i]+binWidth1;
           loc1SegValStart=loc1SegValEnd;
           loc1SegValEnd=val1[i];
        }
        else
        {
	  // Next segment in dhist2 is smaller than the next segment in dhist1
           j+=1;
           loc2SegStart=loc2[j];
           loc2SegEnd=loc2[j]+binWidth2;
           loc2SegValStart=loc2SegValEnd;
           loc2SegValEnd=val2[j];
        }
	// Move the start of the segment under consideration 
        curStartVal=curEndVal;
	// Move the start of the x value under consideration 
        loc2Start=loc2End;
        loc1Start=loc1End;
   }
    return res;
}

