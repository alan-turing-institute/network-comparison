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
//    jumpValues[0]=constantVersion(loc1+1,val1,loc1,val1);
//    jumpValues[1]=constantVersion(loc1+0.5,val1,loc1,val1);
//    jumpValues[2]=constantVersion(loc1+0.25,val1,loc1,val1);
//    jumpValues[3]=constantVersion(loc1+0.125,val1,loc1,val1);
//    jumpValues[4]=constantVersion(loc1+0.06125,val1,loc1,val1);
//    jumpValues[5]=constantVersion(loc1+0.03,val1,loc1,val1);
//    jumpValues[6]=constantVersion(loc1+2,val1,loc1,val1);
//    jumpValues[5]=constantVersion(loc1+0.03,val1,loc1,val1);

    double jumpV1[8];
    jumpV1[0]=2.0;
    jumpV1[1]=1.0;
    jumpV1[2]=0.5;
    jumpV1[3]=0.25;
    jumpV1[4]=0.125;
    jumpV1[5]=0.0625;
    jumpV1[6]=0.015625;
    jumpV1[7]=0.0078125;
    double jumpV2[8];
    for (j=0;j<8;j++)
    {jumpV2[j]=constantVersion(loc1+jumpV1[j],val1,loc1,val1);}
    for (j=0;j<8;j++)
    {
        temp1=constantVersion(loc1+jumpV1[j],val1,loc1,val1);
        if (temp1<jumpV2[j])
        {jumpV2[j]=temp1;}
    }

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
    int count2=0;
    int count3=0;
    double temp2;
    for (i=0;i<offsets.size();i++)
    {
        offset=offsets[i];
        if (offset<offsetLimit)
        {
          //  if (offset<0.01)
           // {std::cout << offset << " skipped " << offsetLimit << "\n";}
            count+=1;
            count2+=1;
            continue;
        }
        offsetDiff=abs(offset-prevOffset);
        // could improve this step by just defining the next valid point
        // but for now lets not do that.
        //
        count1+=1;
        count3+=1;
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
       if (temp1<jumpV2[7])
       {continue;}
       offsetLimit=offset;
       for (j=0;j<8;j++)
       {
           temp2=floor(temp1/jumpV2[j]);
           temp1-=temp2*jumpV2[j];
           offsetLimit+=temp2*jumpV1[j];
       }
      //  std::cout << offset << " " << offsetLimit << " " << count2 << " " << count3 <<" run with\n";
        count2=0;
        count3=0;
    }
    std::cout << " i saved " << count << " calls to the emd function out of " << count1 << "\n";
//    std::cout << " result " << res << " offset " << bestOffset<< "\n";
    return res;
}

