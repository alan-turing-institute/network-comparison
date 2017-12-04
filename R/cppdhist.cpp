// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
using namespace Rcpp;

//double abs(double x)
//{
//    if (x>0)
//    {return x;}
//    else
//    {return -x;}
//}

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
/*
// [[Rcpp::export]]
double constantVersionOnlyCentre(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
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
    std::vector<double> offsets(offsets1.begin(),offsets1.end());
    std::sort(offsets.begin(),offsets.end());
    int i;
    for (i=0;i<offsets.size();i++)
    {

    }
}
*/
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

    double jumpOffsets[8];
    jumpOffsets[0]=2.0;
    jumpOffsets[1]=1.0;
    jumpOffsets[2]=0.5;
    jumpOffsets[3]=0.25;
    jumpOffsets[4]=0.125;
    jumpOffsets[5]=0.0625;
    jumpOffsets[6]=0.015625;
    jumpOffsets[7]=0.0078125;
    double jumpMinSelfEmds[8];
    for (j=0;j<8;j++)
    {
      jumpMinSelfEmds[j]=constantVersion(loc1+jumpOffsets[j],val1,loc1,val1);
    }
    // for (j=0;j<8;j++)
    // {
    //     temp1=constantVersion(loc1+jumpOffsets[j],val1,loc1,val1);
    //     if (temp1<jumpMinSelfEmds[j])
    //     {jumpMinSelfEmds[j]=temp1;}
    // }

    double bestEmd;
    // Initial guess is offset zero. Dhists are min aligned, so this is actually probably
    // a reasonable starting point
    double bestOffset=0;
    bestEmd=constantVersion(loc1,val1,loc2,val2);
    double prevValue=0;
    double offset;
    double offsetLimit=offsets[0]-1;
    double currentEmd;
    double emdDifference;
    int skippedEmdCalls=0;
    int evaluatedEmdCalls=0;
    int count2=0;
    int count3=0;
    double jumpQuotient;
    for (i=0;i<offsets.size();i++)
    {
        offset=offsets[i];
        if (offset<offsetLimit)
        {
          //  if (offset<0.01)
           // {std::cout << offset << " skipped " << offsetLimit << "\n";}
           skippedEmdCalls+=1;
           count2+=1;
            continue;
        }
        // could improve this step by just defining the next valid point
        // but for now lets not do that.
        //
        //skippedEmdCalls+=1;
       evaluatedEmdCalls+=1;
        count3+=1;
        // okay we have to make the expensive call
        //
        currentEmd=constantVersion(loc1+offset,val1,loc2,val2);
        if (currentEmd<bestEmd)
        {
            bestEmd=currentEmd;
            bestOffset=offset;
        }
       // lets work out when the next valid offset will be.
       emdDifference=(currentEmd-bestEmd); // this is the amount we have to deal this.
       if (emdDifference<jumpMinSelfEmds[7])
       {continue;}
       offsetLimit=offset;
       for (j=0;j<8;j++)
       {
           jumpQuotient=floor(emdDifference/jumpMinSelfEmds[j]);
           emdDifference-=jumpQuotient*jumpMinSelfEmds[j];
           offsetLimit+=jumpQuotient*jumpOffsets[j];
       }
        std::cout << offset << " " << offsetLimit << " " << count2 << " " << count3 <<" run with\n";
        count2=0;
        count3=0;
    }
    std::cout << " i saved " << skippedEmdCalls << " calls to the emd function out of " << skippedEmdCalls+evaluatedEmdCalls << " (" << (double)skippedEmdCalls / (double)(skippedEmdCalls+evaluatedEmdCalls) << ")\n";
//    std::cout << " result " << res << " offset " << bestOffset<< "\n";
    return bestEmd;
}

