// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
using namespace Rcpp;


// stores the offset information
// could swap this to be a different data structure
// perhaps a array of vectors
// would then pay a small amount of reallocation
// but would replace while loop - therefore O(num_offsets)
// if statements with a for loop.
struct offsetInfo{
    double offset;
    int loc1loc;
    bool operator < (const offsetInfo& s1) const
    {
        if (offset<s1.offset)
        {
            return 1;
        }
        else if (offset==s1.offset)
        {
            return (loc1loc<s1.loc1loc);
        }
        else
        {
            return 0;
        }
    }
};

// stores the current state of each bin
struct binInfo {
    // starting location
    double loc2Start;
    // end location
    double loc2End;
    // starting height loc2
    double loc2Height;
    // starting height loc1
    double loc1Height;
    // starting index loc1
    int loc1IndexStart;
    // ending index loc1
    int loc1IndexEnd;
};


double updateBin(struct binInfo data,NumericVector val1)
{
    // computes the new emd difference
    // per unit offset
    // for a given bin
    double diff1;
    double h1l1;
    double v1;
    // if there is points in the bin
    // compute the difference between the
    // beginning and end bin
    if (data.loc1IndexEnd!=data.loc1IndexStart)
    {
        h1l1=data.loc1Height;
        v1=data.loc2Height;
        diff1=std::abs(h1l1-v1)-std::abs(val1[data.loc1IndexEnd-1]-v1);
    }
    else
    {diff1=0;}
    return diff1;
}

double fastIterator(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2, std::vector<struct offsetInfo> offsets)
{
    int i,j,k;
    // changes the position of loc1 to place it at the
    // smallest offset
    // Could remove this with the cost of adding the offset to different places.
      loc1=loc1+offsets[0].offset; // this modifies the vector in place this could be a problem

      //The smallest offset
      double minOffset=offsets[0].offset;

      // Update the offset vector with the new offsets
      // again could replace this with a lot of minOffset addition.
      for (i=0;i<offsets.size();i++)
      {
        offsets[i].offset-=minOffset; //update the offsets
      }

      // A vector which store which bin each point
      // on loc1 is in.
      std::vector<int> loc1toBin(loc1.size());


      //A arrays of structs which stores all
      //of the information for each bin
     struct binInfo ls1[loc2.size()+1];

     //initial the locations on the struct
     //Could replace with this a default value.
     for (i=0;i<loc2.size()+1;i++)
     {
         ls1[i].loc1IndexStart=-1;
         ls1[i].loc1IndexEnd=-1;
     }
     // Currently a magic constant should replace with with one that is data dependent
      ls1[0].loc2Start=-10000;
      //the end of the 0th bin
      ls1[0].loc2End=loc2[0];
      //loc2 height at the beginning of the 0th bin
      ls1[0].loc2Height=0;
      //loc1 height at the beginning of the 0th bin
      ls1[0].loc1Height=0;
      // index of the first element in this bin
    ls1[0].loc1IndexStart=0;

//    for i in range(len(loc2)-1):
    for (i=0;i<loc2.size()-1;i++)
    {
//        set the start location of the bin
          ls1[i+1].loc2Start=loc2[i];
//        set the end location of the bin
          ls1[i+1].loc2End=loc2[i+1];
//        set the height of the bin
          ls1[i+1].loc2Height=val2[i];
    }
    // Setup the final bin
      ls1[loc2.size()].loc2Start=loc2[loc2.size()-1];
      // magic constant for end of the bin
      // should be changed to something data dependent
      ls1[loc2.size()].loc2End=1000;
      // set the height of the final bin
      ls1[loc2.size()].loc2Height=val2[val2.size()-1];

    // okay lets add the remaining fields to the data structure.
    j=0;
    ls1[0].loc1IndexStart=0;
    for (i=0;i<loc1.size();i++)
    {
        // Swap to next bin when the offseted bin aligns
        while ((loc2[j]<=loc1[i]) || std::abs(loc2[j]-loc1[i])<0.0000000001)
        {
            //move to the next bin
            ls1[j].loc1IndexEnd=i;
            j+=1;
            if (j==loc2.size())
            {break;}
            ls1[j].loc1IndexStart=i;
            ls1[j].loc1Height=val1[i-1];
        }
        // set the bin for loc1 i
        loc1toBin[i]=j;
        //annoying pattern
        //could do with a goto here!
        //Shouldnt use that much time but would
        //be good to know the better pattern
        if (j==loc2.size())
        {break;}
    }
    // need to deal with the remaining bins
    ls1[j].loc1IndexEnd=i;

    // this is setting the final bin.
    ls1[loc2.size()+1-1].loc1IndexStart=i-1; // this line will be a problem maybe?
    ls1[loc2.size()+1-1].loc1IndexEnd=loc1.size()-1;

    // adding the loc1 locs for each bin
    for (k=i;k<loc1.size()-1;i++)
    {loc1toBin[k]=loc2.size()+1-1;}

    // fill the details of the remaining bins
    for (k=j+1;k<loc2.size()+1;k++)
    {
        ls1[k].loc1IndexStart=loc1.size();
        ls1[k].loc1IndexEnd=loc1.size();
        ls1[k].loc1Height=1.0;
    }

    // fix the first bins height
    // in theory we should be able to
    // remove this?
    ls1[0].loc1Height=0.0;

//    ## first pass on the data.
//
//    Store the current result
    double res=0;

    // vector which stores the diffs for each
    // bin to speed up updates
    std::vector<double> diffs(loc2.size()+1);
    double diff1=0;
    double s1;
    double temp1;
    double se;
    double v1;
    double h1l1;
    double prevH;
    // loop over bins
    // this constructs the first res
    // and constructs the diff vector
    for (i=0;i<loc2.size()+1;i++)
    {
        s1=ls1[i].loc2Start;
        se=ls1[i].loc2End;
        v1=ls1[i].loc2Height;
        h1l1=ls1[i].loc1Height;
        // if there is multiple points in the bin
        if (ls1[i].loc1IndexEnd!=ls1[i].loc1IndexStart)
        {
            prevH=h1l1;
            // iterate over the multiple point
            for (j=ls1[i].loc1IndexStart;j<ls1[i].loc1IndexEnd;j++)
            {

                // add this to the bin
                res+=(loc1[j]-s1)*std::abs(v1-prevH);
                if (res<0)
                {
                    std::cout << "The routine has failed, please submit a bug report" << "\n";
                    return -1;
                }
                // update to the end of the previous bin
                prevH=val1[j];
                s1=loc1[j];
            }
//            # final point
            res+=(se-loc1[ls1[i].loc1IndexEnd-1])*std::abs(val1[ls1[i].loc1IndexEnd-1]-v1); // could be an out by 1 error here.
            diffs[i]=std::abs(h1l1-v1)-std::abs(val1[ls1[i].loc1IndexEnd-1]-v1);
        }
        else
        {
            // case where there no points in the bin
            res+=(double)(ls1[i].loc2End-ls1[i].loc2Start)*(double)std::abs(ls1[i].loc2Height-ls1[i].loc1Height);
            diffs[i]=0;
        }
    }

    //completed first pass, now lets continue
    //
    //
    // Best res and offset
    double curBestEmd=res;
    double curBestOffset=minOffset;
    // variable which holds the current
    // diff per unit offset
    double diffs1=0;
    for (i=0;i<diffs.size();i++)
    {diffs1+=diffs[i];}
    // The old offset
    // need to compute the offset difference
    double oldOffset=0;
    // the bin that
    // needs to be updated
    int bin2Update;
    // increase i
    // to that we start with an offset
    // that is not the first one
    // (we have already computed the first offset)
    for (i=0;i<offsets.size();i++)
    {
        if (offsets[i].offset!=0)
        {break;}
    }
    // iterate over the remaining bins
    // TODO Should this be less than or equal to?
    while (i<offsets.size())
    {
        // update res with the new offset
        res+=(double)(offsets[i].offset-oldOffset)*diffs1;

        // update best res
        // if better
        if (res<curBestEmd)
        {
            curBestEmd=res;
            curBestOffset=offsets[i].offset;
        }

        //update old offset
        oldOffset=offsets[i].offset;

        // Iterate over locations that move bin with this offset
        while (std::abs(offsets[i].offset-oldOffset)<0.000000001) // not massively happy with this tbh
        {
            // The bin that needs to update
            bin2Update=loc1toBin[offsets[i].loc1loc];

            // update the bin that this node is in
            loc1toBin[bin2Update]+=1;
            // move the node (assumed to be the last node)
            // to the next bin
            ls1[bin2Update].loc1IndexEnd-=1;
            ls1[bin2Update+1].loc1IndexStart-=1;
            ls1[bin2Update+1].loc1Height=val1[ls1[bin2Update+1].loc1IndexStart-1]; // i think that this right??

            // update the diffs for the new bin
            diff1=updateBin(ls1[bin2Update],val1);
            diffs1+=diff1-diffs[bin2Update];
            diffs[bin2Update]=diff1;

            //update the diffs for the following bin
            diff1=updateBin(ls1[bin2Update+1],val1);
            diffs1+=diff1-diffs[bin2Update+1];
            diffs[bin2Update+1]=diff1;

            // increase the offset considered
            i+=1;
        }
    }
  loc1=loc1-offsets[0].offset; // this modifies the vector in place this could be a problem
    std::cout << curBestEmd << " " << curBestEmd <<"\n";
    return curBestEmd;
}

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


// [[Rcpp::export]]
double constantVersionExhaustiveCombined(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
    int i,j;
    // this is an experiment if this will be faster than running the standard exhaustive search

    std::vector<struct offsetInfo> offsets;
    //Constructing the offsets
    for (i=0;i<loc1.size();i++)
    {
        for (j=0;j<loc2.size();j++)
        {
            struct offsetInfo product2;
            product2.offset=loc2[j]-loc1[i];
            product2.loc1loc=i;
            offsets.push_back(product2);
        }
    }
    // We need this to sort first by offset and then by the bin id
    // (this is important as it is needed for the updating)
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
    double offsetLimit=offsets[0].offset-1;
    double currentEmd;
    double emdDifference;
    int skippedEmdCalls=0;
    int evaluatedEmdCalls=0;
    int count2=0;
    int count3=0;
    double jumpQuotient;
    i=-1;
    while (i<offsets.size())
    {
        i+=1;i+=1;
        offset=offsets[i].offset;
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
       {
           if (offsets.size()-i>1000)
           {
            // call in the troops
               int numIterations;
                std::vector<struct offsetInfo> newOffsets(offsets.begin()+i,offsets.begin()+i+offsets.size()-1);
                currentEmd=fastIterator(loc1,val1,loc2,val2,newOffsets);
                    i+=offsets.size()-1;
        if (currentEmd<bestEmd)
        {
            bestEmd=currentEmd;
            bestOffset=offset;
        }

           }
           {
               continue;
           }

       }
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

