// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
using namespace Rcpp;


inline double get_segment(double start,double end,double val1_start,double val1_end,double val2_start,double val2_end)
{
  double length;
  length = end - start;
  double topTriangle;
  double topRectangle;
  double bottomTriangle;
  double midPoint;
  double midValue;
  double res=0;
  if (val1_start > val2_start) 
  {
    if (val1_end >= val2_end) 
    {
     // They are in the same order no bowtie
     // top triangle  
//     std::cout << "\n       Path1";
      topTriangle = 0.5*length*(val1_end-val1_start);
      topRectangle = length*(val1_start-val2_start);
      bottomTriangle = 0.5*length*(val2_end-val2_start); 
      return topTriangle+topRectangle-bottomTriangle;
    }
    else
    {
//        std::cout << "\n       Path2";
        //bowtie
        // lets make this really simple as the compiler 
        // will combine the expressions as needed
        midPoint = (val1_start-val2_start)/((val2_end-val2_start) - (val1_end-val1_start));
        midValue = val1_start + midPoint*(val1_end-val1_start);
        midPoint = midPoint*length;
//        std::cout << "\n       midPoint: " << midPoint << " midValue: " << midValue << "\n";
        
        topTriangle = 0.5*midPoint*(midValue-val1_start);
        topRectangle = midPoint*(val1_start-val2_start);
        bottomTriangle = 0.5*midPoint*(midValue-val2_start); 
        
        res = topTriangle+topRectangle-bottomTriangle;
        
        topTriangle = 0.5*(length-midPoint)*(val2_end-midValue);
        topRectangle = 0; // midPoint*(val1_start-val2_start);
        bottomTriangle = 0.5*(length - midPoint)*(val1_end - midValue); 
        res += topTriangle+topRectangle-bottomTriangle;
        return res;
    }
  }
  else
    {
      if (val1_end > val2_end) 
      {
//     std::cout << "\n       Path3";
        //bowtie
        midPoint = (val2_start-val1_start)/((val1_end-val1_start) - (val2_end-val2_start));
        midValue = val2_start + midPoint*(val2_end-val2_start);
        midPoint = midPoint*length;
//        std::cout << "\n       midPoint: " << midPoint << " midValue: " << midValue << "\n";
        
        topTriangle = 0.5*midPoint*(midValue-val2_start);
        topRectangle = midPoint*(val2_start-val1_start);
        bottomTriangle = 0.5*midPoint*(midValue-val1_start); 
        
        res = topTriangle+topRectangle-bottomTriangle;
        
        topTriangle = 0.5*(length-midPoint)*(val1_end-midValue);
        topRectangle = 0; // midPoint*(val1_start-val2_start);
        bottomTriangle = 0.5*(length - midPoint)*(val2_end - midValue); 
        res += topTriangle+topRectangle-bottomTriangle;
        return res;
        
      }
      else // same order
      {
//        std::cout << "\n       Path4";
        topTriangle = 0.5*length*(val2_end-val2_start);
        topRectangle = length*(val2_start-val1_start);
        bottomTriangle = 0.5*length*(val1_end-val1_start); 
        return topTriangle+topRectangle-bottomTriangle;
      }
    }
}

inline double get_segment_constrained(double start,double end, double seg1L1, double seg1L2, double seg2L1, double seg2L2, double seg1V1, double seg1V2, double seg2V1, double seg2V2)  
{
            //We have a valid range
            double valStart1, valEnd1, valStart2, valEnd2;   
            valStart1 = seg1V1 + (seg1V2-seg1V1)*(start - seg1L1)/(seg1L2 - seg1L1);
            valEnd1   = seg1V1 + (seg1V2-seg1V1)*(end   - seg1L1)/(seg1L2 - seg1L1);
            valStart2 = seg2V1 + (seg2V2-seg2V1)*(start - seg2L1)/(seg2L2 - seg2L1);
            valEnd2   = seg2V1 + (seg2V2-seg2V1)*(end   - seg2L1)/(seg2L2 - seg2L1);
            return get_segment(start,end,valStart1,valEnd1,valStart2,valEnd2);
}

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
double NetEmdSmoothV2(NumericVector loc1,NumericVector val1,double binWidth1,NumericVector loc2,NumericVector val2,double binWidth2)
{
  int index1,index2;
  
  double res=0;
  
  // Hist 1
  double curSeg1Loc1; 
  double curSeg1Loc2; 
  double curSeg1Val1; 
  double curSeg1Val2; 
  
  // Hist 1
  double curSeg2Loc1; 
  double curSeg2Loc2; 
  double curSeg2Val1; 
  double curSeg2Val2; 
  
  // need to iterate through regions of constant gradient
  int i123;
  double tempStart;
  double tempEnd;
  double secondStart=-1;
  
  double maxLoc = std::max(loc1[loc1.size()-1] +binWidth1,loc2[loc2.size()-1]+binWidth2 );
  double minLoc = std::min(loc1[0],loc2[0]);
 
  if (loc2[0]<loc1[0]) 
  {
    
    curSeg1Loc1=std::min(loc1[0],loc2[0]); 
    curSeg1Loc2=loc1[0]; 
    curSeg1Val1=0;
    curSeg1Val2=0;
    for (index2=secondStart;index2<loc2.size();index2++) 
    {
              if (index2==-1)
              {
                curSeg2Loc1=minLoc;
                curSeg2Loc2=loc2[0]; 
                curSeg2Val1=0;
                curSeg2Val2=0;
              }
              else
              {
                  curSeg2Loc1=loc2[index2]; 
                  curSeg2Loc2=loc2[index2]+binWidth2; 
                  if (index2==0)
                    {curSeg2Val1=0;}
                  else
                    {curSeg2Val1=val2[index2-1];}
                  curSeg2Val2 = val2[index2]; 
              }
            
            tempStart = std::max(curSeg1Loc1,curSeg2Loc1); 
            tempEnd   = std::min(curSeg1Loc2,curSeg2Loc2); 
            if (tempStart<tempEnd)
            {
              res += get_segment_constrained(tempStart,tempEnd,curSeg1Loc1, curSeg1Loc2, curSeg2Loc1, curSeg2Loc2, curSeg1Val1, curSeg1Val2, curSeg2Val1, curSeg2Val2); 
            }
            if (index2==-1)
              {continue;}
            curSeg2Loc1=loc2[index2]+binWidth2; 
            if (index2==loc2.size()-1)
              {curSeg2Loc2= maxLoc;}
            else
              {curSeg2Loc2=loc2[index2+1];}
            curSeg2Val1=val2[index2]; 
            curSeg2Val2=val2[index2]; 
        
            tempStart = std::max(curSeg1Loc1,curSeg2Loc1); 
            tempEnd   = std::min(curSeg1Loc2,curSeg2Loc2); 
            if (tempStart<tempEnd)
            {
              res += get_segment_constrained(tempStart,tempEnd,curSeg1Loc1, curSeg1Loc2, curSeg2Loc1, curSeg2Loc2, curSeg1Val1, curSeg1Val2, curSeg2Val1, curSeg2Val2); 
            }
            if (curSeg1Loc2<curSeg2Loc1)
              {break;}
        }
    }
    
  
  
  
  
  for (index1=0;index1<loc1.size();index1++) 
  {
      for (i123=0;i123<2;i123++) 
      {
            if (i123==0)
            {
              curSeg1Loc1=loc1[index1]; 
              curSeg1Loc2=loc1[index1]+binWidth1; 
              if (index1==0)
                {curSeg1Val1=0;}
              else
                {curSeg1Val1=val1[index1-1];}
            }
            else
            {
              curSeg1Loc1=loc1[index1]+binWidth1; 
              if (index1==loc1.size()-1)
                {curSeg1Loc2=maxLoc; }
              else
                {curSeg1Loc2=loc1[index1+1];}
              curSeg1Val1=val1[index1]; 
            }
            curSeg1Val2=val1[index1]; 
            if (curSeg1Loc1==curSeg1Loc2)
              {continue;}
    for (index2=secondStart;index2<loc2.size();index2++) 
    {
            if (index2>0)
            {
              if (index2<loc2.size()-2)
              {
                  if (loc2[index2+2]+binWidth2<curSeg1Loc1)
                  {
                    secondStart=index2;
                    continue;
                  }
              }
            }
            if (index2==-1)
            {
              curSeg2Loc1=minLoc;
              curSeg2Loc2=loc2[0]; 
              curSeg2Val1=0;
              curSeg2Val2=0;
            }
            else
            {
                curSeg2Loc1=loc2[index2]; 
                curSeg2Loc2=loc2[index2]+binWidth2; 
                if (index2==0)
                  {curSeg2Val1=0;}
                else
                  {curSeg2Val1=val2[index2-1];}
                curSeg2Val2 = val2[index2]; 
            }
          
          tempStart = std::max(curSeg1Loc1,curSeg2Loc1); 
          tempEnd   = std::min(curSeg1Loc2,curSeg2Loc2); 
          if (tempStart<tempEnd)
          {
            res += get_segment_constrained(tempStart,tempEnd,curSeg1Loc1, curSeg1Loc2, curSeg2Loc1, curSeg2Loc2, curSeg1Val1, curSeg1Val2, curSeg2Val1, curSeg2Val2); 
          }
            if (index2==-1)
              {continue;}
            else
            {
                curSeg2Loc1=loc2[index2]+binWidth2; 
                if (index2==loc2.size()-1)
                  {curSeg2Loc2=maxLoc; }
                else
                  {curSeg2Loc2=loc2[index2+1];}
                curSeg2Val1=val2[index2]; 
                curSeg2Val2=val2[index2]; 
            }
          
          tempStart = std::max(curSeg1Loc1,curSeg2Loc1); 
          tempEnd   = std::min(curSeg1Loc2,curSeg2Loc2); 
          if (tempStart<tempEnd)
          {
            res += get_segment_constrained(tempStart,tempEnd,curSeg1Loc1, curSeg1Loc2, curSeg2Loc1, curSeg2Loc2, curSeg1Val1, curSeg1Val2, curSeg2Val1, curSeg2Val2); 
          }
          if (curSeg1Loc2<curSeg2Loc1)
          {
            break;
          }
      }
    }
  }
  return res;
}
