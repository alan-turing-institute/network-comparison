// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
//#include "src/fastNoSmooth.hpp"
#include <cstdio>
#include <ctime>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

double NetEmdConstantMedianVersion(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
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


bool pairComp(std::pair<double,double> v1,std::pair<double,double> v2)
{
    if (v1.first<v2.first)
    {
        return 1;
    }
    return 0;
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
NumericVector NetEmdExhaustiveMedian(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
  //init
   double res=0;
   double curVal1,curVal2;
   double curPos;
   double temp1;
   int count;
   int i,j,k;
   // structure difference and weight
   std::vector< std::pair<double,double> > offsetsWeights;
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
   double curloc1ValStart=0;
   double curloc1ValEnd=val1[0];
   double curloc2ValStart=0;
   double curloc2ValEnd=val2[0];
   curVal1=0;
   curVal2=0;
   // stores the result
   res=0;
   //TODO be worried about adding lots of small numbers

   // current location on hist 1 and hist 2
   i=0;
   j=0;
   double value;
   double offset;
   std::clock_t start;
   double duration;
   start = std::clock();
   while (1)
    {
        if (curloc1ValEnd<curloc2ValEnd)
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc1ValEnd-curloc2ValStart;}
            else
            { value=curloc1ValEnd-curloc1ValStart;}
            offset=loc2[j]-loc1[i];
            if (value>0.0000000000001) // need to discuss this constant
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                offsetsWeights.push_back(temp1);
            }
            i+=1;
            if (i==loc1.size())
            {break;}
            curloc1ValStart=curloc1ValEnd;
            curloc1ValEnd=val1[i];
        }
        else
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc2ValEnd-curloc2ValStart;}
            else
            {value=curloc2ValEnd-curloc1ValStart;}
            offset=loc2[j]-loc1[i];
            if (value>0.0000000000001)
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                offsetsWeights.push_back(temp1);
            }
            j+=1;
            if (j==loc2.size())
            {break;}
            curloc2ValStart=curloc2ValEnd;
            curloc2ValEnd=val2[j];
        }
    }
    if (i<loc1.size())
    {
        for (k=i+1;k<loc1.size();k++)
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc1ValEnd-curloc2ValStart;}
            else
            { value=curloc1ValEnd-curloc1ValStart;}
            offset=loc2[loc2.size()-1]-loc1[k];
            if (value>0.0000000000001)
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                offsetsWeights.push_back(temp1);
            }
            curloc1ValStart=curloc1ValEnd;
            curloc1ValEnd=val1[k];
        }
    }
    else
    {
        for (k=j+1;k<loc2.size();k++)
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc2ValEnd-curloc2ValStart;}
            else
            {value=curloc2ValEnd-curloc1ValStart;}
            offset=loc2[k]-loc1[loc1.size()-1];
            if (value>0.0000000000001)
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                offsetsWeights.push_back(temp1);
            }
            curloc2ValStart=curloc2ValEnd;
            curloc2ValEnd=val2[k];
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"largeLoop: "<< duration <<'\n';
    NumericVector result(2);

    start = std::clock();

    std::sort(offsetsWeights.begin(),offsetsWeights.end(),pairComp);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<< " sort: "<< duration <<'\n';

    start = std::clock();

    double currentWeight=0;
    double curBest;
    count=0;
    NumericVector loc0;
    for (i=0;i<offsetsWeights.size();i++)
    {
    //    loc0=loc1+offsetsWeights[i].first;
    //    curBest=NetEmdConstant(loc0,val1,loc2,val2);
//        std::cout << offsetsWeights[i].first << "," << offsetsWeights[i].second << ","<< currentWeight << "," << curBest <<  "\n";
        value=offsetsWeights[i].second;
        if (currentWeight+value>0.5)
        {
            //this is the point
            currentWeight+=value;
            loc0=loc1+offsetsWeights[i].first;
            NumericVector result(2);
            duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
            std::cout<<"find location:  "<< duration <<'\n';
	    start = std::clock();
            result(0)=NetEmdConstantMedianVersion(loc0,val1,loc2,val2);
            result(1)=offsetsWeights[i].first;
            duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

            std::cout<<"emd call: "<< duration <<'\n';
            return result;
        }
        else
        { currentWeight+=value; }
    }
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
NumericVector NetEmdExhaustiveMedianMk2(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2,double lowerLimit,double upperLimit)
{
  //init
   double res=0;
   double curVal1,curVal2;
   double curPos;
   double temp1;
   int count;
   int i,j,k;
   double lowerSum=0;
   double upperSum=0;
   double lowerSumOutside=0;
   double upperSumOutside=0;
   // structure difference and weight
   std::vector< std::pair<double,double> > offsetsWeightsLower;
   std::vector< std::pair<double,double> > offsetsWeightsUpper;
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
   double curloc1ValStart=0;
   double curloc1ValEnd=val1[0];
   double curloc2ValStart=0;
   double curloc2ValEnd=val2[0];
   curVal1=0;
   curVal2=0;
   // stores the result
   res=0;
   //TODO be worried about adding lots of small numbers

   // current location on hist 1 and hist 2
   i=0;
   j=0;
   double value;
   double offset;
   std::clock_t start;
   double duration;
   start = std::clock();
   int numberOfSegments=0;
   while (1)
    {
        if (curloc1ValEnd<curloc2ValEnd)
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc1ValEnd-curloc2ValStart;}
            else
            { value=curloc1ValEnd-curloc1ValStart;}
            offset=loc2[j]-loc1[i];
            if (value>0.0000000000001) // need to discuss this constant
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                if (offset<0)
                {
		    if (lowerLimit<offset)
		    {
			    offsetsWeightsLower.push_back(temp1);
		    }
                    lowerSum+=value;
		    numberOfSegments+=1;
                }
                else
                {
		    if (upperLimit>offset)
		    {
			    offsetsWeightsUpper.push_back(temp1);
		    }
                    upperSum+=value;
		    numberOfSegments+=1;
                }
            }
            i+=1;
            if (i==loc1.size())
            {break;}
            curloc1ValStart=curloc1ValEnd;
            curloc1ValEnd=val1[i];
        }
        else
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc2ValEnd-curloc2ValStart;}
            else
            {value=curloc2ValEnd-curloc1ValStart;}
            offset=loc2[j]-loc1[i];
            if (value>0.0000000000001)
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                if (offset<0)
                {
		    if (lowerLimit<offset)
		    {
			    offsetsWeightsLower.push_back(temp1);
		    }
                    lowerSum+=value;
		    numberOfSegments+=1;
                }
                else
                {
		    if (upperLimit>offset)
		    {
			    offsetsWeightsUpper.push_back(temp1);
		    }
                    upperSum+=value;
		    numberOfSegments+=1;
                }
            }
            j+=1;
            if (j==loc2.size())
            {break;}
            curloc2ValStart=curloc2ValEnd;
            curloc2ValEnd=val2[j];
        }
    }
    if (i<loc1.size())
    {
        for (k=i+1;k<loc1.size();k++)
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc1ValEnd-curloc2ValStart;}
            else
            { value=curloc1ValEnd-curloc1ValStart;}
            offset=loc2[loc2.size()-1]-loc1[k];
            if (value>0.0000000000001)
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                if (offset<0)
                {
		    if (lowerLimit<offset)
		    {
			    offsetsWeightsLower.push_back(temp1);
		    }
                    lowerSum+=value;
		    numberOfSegments+=1;
                }
                else
                {
		    if (upperLimit>offset)
		    {
			    offsetsWeightsUpper.push_back(temp1);
		    }
                    upperSum+=value;
		    numberOfSegments+=1;
                }
            }
            curloc1ValStart=curloc1ValEnd;
            curloc1ValEnd=val1[k];
        }
    }
    else
    {
        for (k=j+1;k<loc2.size();k++)
        {
            if (curloc1ValStart<curloc2ValStart)
            {value=curloc2ValEnd-curloc2ValStart;}
            else
            {value=curloc2ValEnd-curloc1ValStart;}
            offset=loc2[k]-loc1[loc1.size()-1];
            if (value>0.0000000000001)
            {
                std::pair<double,double> temp1;
                temp1.first=offset;
                temp1.second=value;
                if (offset<0)
                {
		    if (lowerLimit<offset)
		    {
			    offsetsWeightsLower.push_back(temp1);
		    }
                    lowerSum+=value;
		    numberOfSegments+=1;
                }
                else
                {
		    if (upperLimit>offset)
		    {
			    offsetsWeightsUpper.push_back(temp1);
		    }
                    upperSum+=value;
		    numberOfSegments+=1;
                }
            }
            curloc2ValStart=curloc2ValEnd;
            curloc2ValEnd=val2[k];
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"largeLoop: "<< duration <<'\n';
    NumericVector result(2);

    start = std::clock();

double currentWeight;
    double curBest;
    NumericVector loc0;
    if (lowerSum>upperSum)
    {
        std::sort(offsetsWeightsLower.begin(),offsetsWeightsLower.end(),pairComp);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"sorting: "<< duration <<'\n';
        std::cout<<"filter: "<< offsetsWeightsLower.size() << " " << numberOfSegments <<'\n';
        currentWeight=lowerSum;
        start = std::clock();
        for (i=offsetsWeightsLower.size()-1;i>=0;i--)
        {
            value=offsetsWeightsLower[i].second;
            currentWeight-=value;
            if (currentWeight<0.5)
            {
                loc0=loc1+offsetsWeightsLower[i].first;
                NumericVector result(2);
                duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		    std::cout<<"searching: "<< duration <<'\n';
		start = std::clock();
                result(0)=NetEmdConstantMedianVersion(loc0,val1,loc2,val2);
                duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		    std::cout<<"emdComputation: "<< duration <<'\n';
                result(1)=offsetsWeightsLower[i].first;
                return result;
            }
        }
    }
    else
    {
        std::sort(offsetsWeightsUpper.begin(),offsetsWeightsUpper.end(),pairComp);
        std::cout<<"filter: "<< offsetsWeightsUpper.size() << " " << numberOfSegments <<'\n';
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"sorting: "<< duration <<'\n';
        currentWeight=upperSum;
        start = std::clock();
        for (i=0;i<offsetsWeightsUpper.size();i++)
        {
            value=offsetsWeightsUpper[i].second;
            currentWeight-=value;
            if (currentWeight<0.5)
            {
                loc0=loc1+offsetsWeightsUpper[i].first;
                NumericVector result(2);
		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		std::cout<<"searching: "<< duration <<'\n';
		start = std::clock();
                result(0)=NetEmdConstantMedianVersion(loc0,val1,loc2,val2);
                duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		std::cout<<"EmdCall: "<< duration <<'\n';
                result(1)=offsetsWeightsUpper[i].first;
                return result;
            }
        }
    }
}

