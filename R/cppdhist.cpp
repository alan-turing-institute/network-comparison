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


// [[Rcpp::export]]
double constantVersionExhaustive(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
    int i,j;
    // this is an experiment if this will be faster than running the standard exhaustive search
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
    double jumpValues[7];
    double temp1;
    jumpValues[0]=constantVersion(loc1+1,val1,loc1,val1);
    jumpValues[1]=constantVersion(loc1+0.5,val1,loc1,val1);
    jumpValues[2]=constantVersion(loc1+0.25,val1,loc1,val1);
    jumpValues[3]=constantVersion(loc1+0.125,val1,loc1,val1);
    jumpValues[4]=constantVersion(loc1+0.06125,val1,loc1,val1);
    jumpValues[5]=constantVersion(loc1+0.03,val1,loc1,val1);
    jumpValues[6]=constantVersion(loc1+2,val1,loc1,val1);
    jumpValues[5]=constantVersion(loc1+0.03,val1,loc1,val1);


    double res;
    double bestOffset=0;
    res=constantVersion(loc1,val1,loc2,val2);
    double prevValue=0;
    double prevOffset=-1000;
    double offset;
    double offsetLimit=-100000;
    double offsetDiff;
    double tempVal;
    int count=0;
    int count1=0;
    double temp2;
    for (i=0;i<offsets.size();i++)
    {
        offset=offsets[i];
        if (offset<offsetLimit)
        {
          //  if (offset<0.01)
           // {std::cout << offset << " skipped " << offsetLimit << "\n";}
            count+=1;
            continue;
        }
            if ((offset)<0.01)
            {
            }
        offsetDiff=abs(offset-prevOffset);
        // could improve this step by just defining the next valid point
        // but for now lets not do that.
        //
        count1+=1;
        // okay we have to make the expensive call
        //
        tempVal=constantVersion(loc1+offset,val1,loc2,val2);
        if (tempVal<res)
        {
            res=tempVal;
            bestOffset=offset;
        }
        prevOffset=offset;
        prevValue=tempVal;
// lets work out when the next valid offset will be.
       temp1=(prevValue-res); // this is the amount we have to deal this.
       if (temp1<jumpValues[5])
       {continue;}
       temp2=floor(temp1/jumpValues[6]);
       offsetLimit=offset+temp2*2;
       temp1-=temp2*jumpValues[6];
       temp2=floor(temp1/jumpValues[0]);
       offsetLimit+=temp2;
       temp1-=temp2*jumpValues[0];
       temp2=floor(temp1/jumpValues[1]);
       offsetLimit+=temp2*0.5;
       temp1-=temp2*jumpValues[1];
       temp2=floor(temp1/jumpValues[2]);
       offsetLimit+=temp2*0.25;
       temp1-=temp2*jumpValues[2];
       temp2=floor(temp1/jumpValues[3]);
       offsetLimit+=temp2*0.125;
       temp1-=temp2*jumpValues[3];
       temp2=floor(temp1/jumpValues[4]);
       offsetLimit+=temp2*0.06125;
       temp1-=temp2*jumpValues[4];
       temp2=floor(temp1/jumpValues[5]);
       offsetLimit+=temp2*0.03;
       temp1-=temp2*jumpValues[5];
        std::cout << offset << " " << offsetLimit << " run with\n";
    }
    std::cout << " i saved " << count << " calls to the emd function out of " << count1 << "\n";
//    std::cout << " result " << res << " offset " << bestOffset<< "\n";
    return res;
}

