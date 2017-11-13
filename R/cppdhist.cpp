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
double constantVersion(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
   double res=0;
   double curVal1,curVal2;
   double curPos;
   int i,j,k;
   int maxLen;
   maxLen=std::max(loc1.size(),loc2.size());
   if (loc1[0]<loc2[0])
   {
       res+=(loc2[0]-loc1[0])*(val1[0]);
       curPos=loc2[0];
       curVal1=val1[0];
       curVal2=val2[0];
   }
   else
   {
       res+=(loc1[0]-loc2[0])*(val2[0]);
       curPos=loc1[0];
       curVal1=val1[0];
       curVal2=val2[0];
   }
   i=1;
   j=1;
   for (k=0;k<=maxLen;k++)
    {
        if (i==loc1.size())
        {break;}
        if (j==loc2.size())
        {break;}
        if (loc1[i]<loc2[j])
        {
            res+=(loc1[i]-curPos)*abs(curVal1-curVal2);
            curVal1=val1[i];
            curPos=loc1[i];
            i+=1;
        }
        else
        {
            res+=(loc2[j]-curPos)*abs(curVal1-curVal2);
            curVal2=val2[j];
            curPos=loc2[j];
            j+=1;
        }
    }
    if (i<loc1.size())
    {
        for (k=i;k<loc1.size();k++)
        {
            res+=(loc1[k]-curPos)*(1-curVal1);
            curVal1=val1[k];
            curPos=loc1[k];
        }
    }
    else
    {
        for (k=j;k<loc2.size();k++)
        {
            res+=(loc2[k]-curPos)*(1-curVal2);
            curVal2=val2[k];
            curPos=loc2[k];
        }
    }
    return res;
}


//#double func(loc1,val1,loc2,val2)
//#//We can assume that the locations are sorted
//#//therefore we just need to merge them
//#   std:vector<double> x;
//#   std:vector<double> val1;
//#   std:vector<double> val2;
//#   int i,j,k;
//#   int maxLen;
//#   maxLen=max(loc1.size(),loc2.size())
//#    i=0;
//#    j=0;
//#   for (k=0;k<maxLen;k++)
//#    {
//#        if (i==loc1.size())
//#        {break;}
//#        if (j==loc2.size())
//#        {break;}
//#        if (loc1[i]<loc2[j])
//#        {
//#            x.push_back(loc1[i]);
//#            i+=1;
//#        }
//#        else
//#        {
//#            x.push_back(loc1[j]);
//#            j+=1;
//#        }
//#    }
//#    if (i<loc1.size())
//#    {
//#        for (k=i;k<loc1.size();k++)
//#        {
//#            x.push_back(loc1[k]);
//#        }
//#    }
//#    else
//#    {
//#        for (k=j;k<loc2.size();k++)
//#        {
//#            x.push_back(loc2[k]);
//#        }
//#    }
//#
