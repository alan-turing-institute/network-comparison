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

// [[Rcpp::export]]
double constantVersionWithUpdates(NumericVector loc1,NumericVector val1,NumericVector loc2,NumericVector val2)
{
    int i,j,k;
    NumericVector result;
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

            //update bin location
            loc1toBin[offsets[i].loc1loc]+=1;
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
    std::cout << curBestEmd << " " << curBestEmd <<"\n";
    return curBestEmd;
}
