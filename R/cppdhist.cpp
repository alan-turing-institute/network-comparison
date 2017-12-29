// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <ctime>

using namespace Rcpp;

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
    double bestEmd;
    // Initial guess is offset zero. Dhists are min aligned, so this is actually probably
    // a reasonable starting point
    double bestOffset=0;
    bestEmd=constantVersion(loc1,val1,loc2,val2);
    double prevValue=0;
    double offset;
    double offsetLimit=-bestEmd; //offsets[0]-1;
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
//       if (emdDifference<jumpMinSelfEmds[7])
//       {continue;}
       offsetLimit=offset;
       offsetLimit=offset+emdDifference;
//       for (j=0;j<8;j++)
//       {
//           jumpQuotient=floor(emdDifference/jumpMinSelfEmds[j]);
//           emdDifference-=jumpQuotient*jumpMinSelfEmds[j];
//           offsetLimit+=jumpQuotient*jumpOffsets[j];
//       }
//        std::cout << offset << " " << offsetLimit << " " << count2 << " " << count3 <<" run with\n";
        count2=0;
        count3=0;
    }
    //std::cout << " i saved " << skippedEmdCalls << " calls to the emd function out of " << skippedEmdCalls+evaluatedEmdCalls << " (" << (double)skippedEmdCalls / (double)(skippedEmdCalls+evaluatedEmdCalls) << ")\n";
    std::cout << " i saved " << evaluatedEmdCalls  << "\n";
//    std::cout << " result " << res << " offset " << bestOffset<< "\n";
    return bestEmd;
}
















struct Interval
{
    double minValue;
    double start;
    double end;
    double startVal;
    double endVal;
    int startOffset;
    int endOffset;
    bool operator < (const Interval& s1) const
    {
        if (minValue<s1.minValue)
        {
            return 1;
        }
        else if (minValue==s1.minValue)
        {
            return (start<s1.start);
        }
        else
        {
            return 0;
        }
    }
};


struct emdDualResult
{
    double leftAlignment;
    double rightAlignment;
    double leftValue;
    double rightValue;
};


//NumericVector emdDualAlignments(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)


struct emdDualResult emdDualAlignments(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
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
   curVal1=0;
   curVal2=0;
   // stores the result
   res=0;

   // current location on hist 1 and hist 2
   i=0;
   j=0;
    double diffs=0;
    double minDiffLeft=1000;
    double minDiffRight=1000;
    double tempDouble;
    // deal with the start diff.
    while (1)
    {
    //        std::cout << "i=" << i;
    //        std::cout << " j=" << j;
    //        std::cout << " res=" << res;
    //        std::cout << "\n";
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
            // so we have a j update so we need to compute the differences.
            if (i>1)
            {
                tempDouble=loc2[j]-loc1[i-1];
                if (tempDouble<minDiffLeft)
                {
                    minDiffLeft=tempDouble;
                }
            }
            tempDouble=loc1[i]-loc2[j];
            if (tempDouble<minDiffRight)
            {
                minDiffRight=tempDouble;
            }
            temp1=(loc2[j]-curPos)*std::abs(curVal1-curVal2);
            res+=temp1;
            diffs+=std::abs(curVal1-val2[j])-std::abs(curVal1-curVal2);
//            std::cout << "i=" << i;
//            std::cout << " j=" << j;
//            std::cout << " diff=" << std::abs(curVal1-val2[j])-std::abs(curVal1-curVal2);
//            std::cout << " diffs=" << diffs;
//            std::cout << " curVal1=" << curVal1;
//            std::cout << " val2[j]=" << val2[j];
//            std::cout << " curVal2=" << curVal2;
//            std::cout << "\n";
            curVal2=val2[j];
            curPos=loc2[j];
 //p           << "\tj=" << j << "\tdiffs=" << diffs <<  "\t " << std::abs(curVal1-val2[j+1])-std::abs(curVal1-val2[j]) << "\n";
            j+=1;
        }
    }
//    std::cout << "res=" << res << " \n\n";
//    std::cout << "gap j=" << j << "\n\n";
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
            diffs+=std::abs(curVal1-val2[k])-std::abs(curVal1-curVal2);
//            std::cout << "i=" << i;
//            std::cout << " j=" << k;
//            std::cout << " diff=" << std::abs(curVal1-val2[j])-std::abs(curVal1-curVal2);
//            std::cout << " diffs=" << diffs;
//            std::cout << " curVal1=" << curVal1;
//            std::cout << " val2[j]=" << val2[k];
//            std::cout << " curVal2=" << curVal2;
//            std::cout << "\n";
            curVal2=val2[k];
            curPos=loc2[k];
 //           std::cout << "i=" << i << "\tj=" << k << "\tdiffs=" << diffs <<  "\t " << std::abs(curVal1-val2[k+1])-std::abs(curVal1-val2[k]) << "\n";
        }
    }
//    std::cout << "res=" << res << " \n\n";
//    NumericVector xx(5);
//    xx(0)=diffs;
//    xx(1)=minDiffLeft;
//    xx(2)=minDiffRight;
//    xx(3)=res+diffs*minDiffLeft;
//    xx(4)=res-diffs*minDiffRight;
//    return xx;
    struct emdDualResult returnObject;
    returnObject.leftAlignment=minDiffLeft;
    returnObject.rightAlignment=minDiffRight;
    returnObject.leftValue=res+diffs*minDiffLeft;
    returnObject.rightValue=res-diffs*minDiffRight;
    return returnObject;
}



// [[Rcpp::export]]
double constantVersionExhaustiveHalf(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
//    ProfilerStart("./eigen-prof.log");
    std::time_t startTime = std::time(nullptr);

    int i,j;
    // this is an experiment if this will be faster than running the standard exhaustive search
    std::unordered_set<double> offsets1;
    int pos=0;

    for (i=0;i<loc2.size();i++)
    {
        for (j=0;j<loc1.size();j++)
        {
            offsets1.insert(loc2[i]-loc1[j]);
        }
    }
    // probably is better way to do this.
    std::time_t setTime = std::time(nullptr);
    std::vector<double> offsets(offsets1.begin(),offsets1.end());
    std::time_t vectTime = std::time(nullptr);
    std::sort(offsets.begin(),offsets.end());



    std::time_t mediumTime = std::time(nullptr);
    std::priority_queue<Interval> priorityQ1;

    double jumpValues[7];
    double temp1;

    struct Interval inter1;
    inter1.start=loc2[1]-loc1[loc1.size()];
    inter1.end=loc2[loc2.size()]-loc1[1];
    inter1.startVal=loc1[loc1.size()]-loc2[1];
    inter1.endVal=loc2[loc2.size()]-loc1[1];
    inter1.minValue=0;
    inter1.startOffset=0;
    inter1.endOffset=offsets.size()-1;
    priorityQ1.push(inter1);

    double prevValue=0;
    double offset;
    double currentEmd;
    double emdDifference;
    int skippedEmdCalls=0;
    int evaluatedEmdCalls=0;
    int count2=0;
    int count3=0;
    double jumpQuotient;
    struct Interval currentItem;
    double midPoint;
    double leftPoint,leftValue,rightPoint,rightValue;
    double t1;
    double midIndex;
    double curBest=loc1[loc1.size()]-loc2[1];
    int count=0;
    int emdCalls=0;
    while (priorityQ1.size()>0)
    {
        count+=1;
        currentItem=priorityQ1.top();
        priorityQ1.pop();
        if (-currentItem.minValue>curBest)
        {
            break;
        }
//        if (count==100)
//        {break;}

        if (currentItem.endOffset-currentItem.startOffset<5)
        {
            for (i=currentItem.startOffset+1;i<currentItem.endOffset;i++)
            {
                currentEmd=constantVersion(loc1+offsets[i],val1,loc2,val2);
                emdCalls+=1;
                if (currentEmd < curBest)
                {
                    curBest=currentEmd;
                }
            }
            continue;
        }
        midIndex=currentItem.startOffset+(currentItem.endOffset-currentItem.startOffset)/2;
        midPoint=offsets[midIndex];
        currentEmd=constantVersion(loc1+midPoint,val1,loc2,val2);
        emdCalls+=1;
        if (currentEmd < curBest)
        {
            curBest=currentEmd;
        }


        leftPoint=currentItem.start;
        rightPoint=midPoint;
        leftValue=currentItem.startVal;
        rightValue=currentEmd;
        t1=(leftValue+rightValue)/2.0+(leftPoint-rightPoint)/2.0;
        struct Interval newInterval1;
        newInterval1.minValue=-t1;
        newInterval1.start=leftPoint;
        newInterval1.end=rightPoint;
        newInterval1.startVal=leftValue;
        newInterval1.endVal=rightValue;
        newInterval1.startOffset=currentItem.startOffset;
        newInterval1.endOffset=midIndex;
        if (newInterval1.start!=newInterval1.end)
        {
        priorityQ1.push(newInterval1);
        }
//        std::cout << "t1=" << t1 << "\n";

        leftPoint=midPoint;
        rightPoint=currentItem.end;
        leftValue=currentEmd;
        rightValue=currentItem.endVal;
        t1=(leftValue+rightValue)/2.0+(leftPoint-rightPoint)/2.0;
        struct Interval newInterval2;
        newInterval2.minValue=-t1;
        newInterval2.start=leftPoint;
        newInterval2.end=rightPoint;
        newInterval2.startVal=leftValue;
        newInterval2.endVal=rightValue;
        newInterval2.startOffset=midIndex;
        newInterval2.endOffset=currentItem.endOffset;
 //       std::cout << leftPoint << " "<< rightPoint << "t1=" << t1 << "\n";
        if (newInterval2.start!=newInterval2.end)
        {
        priorityQ1.push(newInterval2);
        }
    }
    std::cout << "i made " << emdCalls << "\n";
//    std::cout << " result " << res << " offset " << bestOffset<< "\n";
//        ProfilerStop();
//
    std::cout << "set step: " << setTime-startTime << "\n";
    std::cout << "vec step: " << vectTime-setTime << "\n";
    std::cout << "sort step: " << mediumTime-vectTime<< "\n";
    std::cout << "remaining time: " << std::time(nullptr) - mediumTime << "\n";

    return curBest;
}


struct IntervalNoOffset
{
    double minValue;
    double start;
    double end;
    double startVal;
    double endVal;
    bool operator < (const IntervalNoOffset& s1) const
    {
        if (minValue<s1.minValue)
        {
            return 1;
        }
        else if (minValue==s1.minValue)
        {
            return (start<s1.start);
        }
        else
        {
            return 0;
        }
    }
};

// [[Rcpp::export]]
double constantVersionExhaustiveHalfNoOffsetCalc(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
    //std::unordered_set<double> offsets1;
    //int i,j;
    //for (i=0;i<loc2.size();i++)
    //{
    //    for (j=0;j<loc1.size();j++)
    //    {
    //        offsets1.insert(loc2[i]-loc1[j]);
    //    }
    //}
    //std::vector<double> offsets(offsets1.begin(),offsets1.end());
    //std::sort(offsets.begin(),offsets.end());
    std::cout << "im starting\n";
    std::time_t startTime = std::time(nullptr);

    // probably is better way to do this.
    std::time_t setTime = std::time(nullptr);
    std::time_t vectTime = std::time(nullptr);


    std::time_t mediumTime = std::time(nullptr);
    std::priority_queue<IntervalNoOffset> priorityQ1;

    double jumpValues[7];
    double temp1;

    struct IntervalNoOffset inter1;
    inter1.start=loc2[0]-loc1[loc1.size()-1];
    inter1.end=loc2[loc2.size()-1]-loc1[0];
    inter1.startVal=loc1[loc1.size()-1]-loc2[0];
    inter1.endVal=loc2[loc2.size()-1]-loc1[0];
    inter1.minValue=0;
    priorityQ1.push(inter1);

    double prevValue=0;
    double offset;
    double currentEmd;
    double emdDifference;
    int skippedEmdCalls=0;
    int evaluatedEmdCalls=0;
    int count2=0;
    int count3=0;
    double jumpQuotient;
    struct IntervalNoOffset currentItem;
    double midPoint;
    double leftPoint,leftValue,rightPoint,rightValue;
    double t1;
    double midIndex;
    double curBest=loc1[loc1.size()-1]-loc2[0];
    int count=0;
    int emdCalls=0;
    struct emdDualResult emdDualStructureStore;
    while (priorityQ1.size()>0)
    {
        count+=1;
        currentItem=priorityQ1.top();
        priorityQ1.pop();
        if (-currentItem.minValue>curBest)
        {
            std::cout << "currentItemMinValue=" << currentItem.minValue;
            std::cout << " curBest=" << curBest;
            std::cout << "\ni hit the break\n";
            break;
        }
//        if (count==100)
//        {break;}
        midPoint=currentItem.start+(currentItem.end-currentItem.start)/2;
//        std::cout << "midPoint=";
//        std::cout << midPoint;
//        std::cout << " start=";
//        std::cout << currentItem.start;
//        std::cout << " end=";
//        std::cout << currentItem.end;
//        std::cout << "\n";
        emdDualStructureStore=emdDualAlignments(loc1+midPoint,val1,loc2,val2);
    //    double leftClose=1000;
    //    double rightClose=1000;
    //    for (i=0;i<offsets.size();i++)
    //    {
    //        if (offsets[i]<midPoint)
    //        {
    //            leftClose=offsets[i];
    //        }
    //        else
    //        {
    //            rightClose=offsets[i];
    //            break;
    //        }
    //    }
    //    std::cout << "leftClose=";
    //    std::cout << leftClose;
    //    std::cout << "rightClose=";
    //    std::cout << rightClose;
    //    std::cout << "\n";
    //    std::cout << " leftValue=";
    //    std::cout << emdDualStructureStore.leftValue;
    //    std::cout << " ";
    //    std::cout << " rightValue=";
    //    std::cout << emdDualStructureStore.rightValue;
    //    std::cout << " ";
    //    std::cout << " leftAlignment=";
    //    std::cout << emdDualStructureStore.leftAlignment;
    //    std::cout << " ";
    //    std::cout << " rightAlignment=";
    //    std::cout << emdDualStructureStore.rightAlignment;
    //    std::cout << "\n";
    //    std::cout << constantVersion(loc1+midPoint+emdDualStructureStore.leftAlignment,val1,loc2,val2);
    //    std::cout << "\n";
    //    std::cout << constantVersion(loc1+midPoint-emdDualStructureStore.rightAlignment,val1,loc2,val2);
    //    std::cout << "\n";
    ////    if (emdCalls==2)
    //    {

    //    return emdDualStructureStore.leftValue;
    //    }
        emdCalls+=1;
        if (emdDualStructureStore.rightValue < curBest)
        {
            curBest=emdDualStructureStore.rightValue;
        }
        if (emdDualStructureStore.leftValue < curBest)
        {
            curBest=emdDualStructureStore.leftValue;
        }


        leftPoint=currentItem.start;
        rightPoint=midPoint-emdDualStructureStore.rightAlignment;
        if ((rightPoint-leftPoint)>0.000000000001)
        {
            leftValue=currentItem.startVal;
            rightValue=emdDualStructureStore.rightValue;
            t1=(leftValue+rightValue)/2.0+(leftPoint-rightPoint)/2.0;
            if (t1<curBest)
            {
                struct IntervalNoOffset newInterval1;
                newInterval1.minValue=-t1;
                newInterval1.start=leftPoint;
                newInterval1.end=rightPoint;
                newInterval1.startVal=leftValue;
                newInterval1.endVal=rightValue;
                priorityQ1.push(newInterval1);
            }
        }
//        std::cout << "t1=" << t1 << "\n";

        leftPoint=midPoint+emdDualStructureStore.leftAlignment;
        rightPoint=currentItem.end;
        leftValue=emdDualStructureStore.leftValue;
        rightValue=currentItem.endVal;
        if ((rightPoint-leftPoint)>0.000000000001)
        {
            t1=(leftValue+rightValue)/2.0+(leftPoint-rightPoint)/2.0;
            if (t1<curBest)
            {
                struct IntervalNoOffset newInterval2;
                newInterval2.minValue=-t1;
                newInterval2.start=leftPoint;
                newInterval2.end=rightPoint;
                newInterval2.startVal=leftValue;
                newInterval2.endVal=rightValue;
                priorityQ1.push(newInterval2);
            }
        }
    }
    std::cout << "i made " << emdCalls << "\n";
//    std::cout << " result " << res << " offset " << bestOffset<< "\n";
//        ProfilerStop();
//
    std::cout << "set step: " << setTime-startTime << "\n";
    std::cout << "vec step: " << vectTime-setTime << "\n";
    std::cout << "sort step: " << mediumTime-vectTime<< "\n";
    std::cout << "remaining time: " << std::time(nullptr) - mediumTime << "\n";

    return curBest;
}
