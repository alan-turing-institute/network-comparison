#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

struct offsetInfo{
    double offset;
    int loc1loc;
}

struct binInfo {
    int loc2Start;
    int loc2End;
    double loc2Height;
    double loc1Height;
    int loc1IndexStart;
    int loc2IndexStart;
}

//def constantVersionWithAccumPythonVersion(loc1,val1,loc2,val2,offsets):
NumericVector constantVersionWithUpdates(loc1,val1,loc2,val2)
{
    int i,j;
    NumericVector result;
    std::vector offsets<struct offsetInfo>;
//    for i in range(len(loc1)):
//        for j in range(len(loc2)):
//            offsets.append([loc2[j]-loc1[i],i,j])
    for (i=0;i<loc1.size();i++)
    {
        for (j=0;i<loc2.size();j++)
        {
            struct offsetInfo product2;
            product2.offset=loc2[j]-loc1[i];
            product2.loc1loc=i;
            offsets.insert(product2);
        }
    }
//    offsets.sort()
    std::sort(offsets.begin(),offsets.end()); // We need this to sort first by offset and then by the bin id (this is important as it is needed for the ordering)
//    for i in range(len(loc1)):
//        loc1[i]+=offsets[0][0]
      loc1+=offsets[0][0]; // this modifies the vector in place this could be a problem
      double minOffset=offsets[0][0];
//    minOffset=offsets[0][0]
//    for item in offsets:
//        item[0]-=temp1
      for (i=0;i<offsets.size();i++)
      {
        offsets[i][0]-=minOffset; //update the offsets
      }
//    loc1toBin={}
      std::vector<int> loc1toBin(loc1.size());
//    ls1=[]
     struct binInfo ls1[loc2.size()+1];
//
//    t1['loc2start']=-10000
      ls1[0].loc2Start=-10000;
//    t1['loc2end']=loc2[0]
      ls1[0].loc2End=loc2[0];
//    t1['loc2Height']=0
      ls1[0].loc2Height=0;
//    t1['loc1Height']=getHeightAtPoint(loc1,val1,-10000)
      ls1[0].loc1Height=0;
//    t1['otherPoints'],temp1=getOtherPoints(loc1,val1,-10000,loc2[0])
    ls1[0].loc1IndexStart=0
//    for j in temp1:
//        loc1toBin[j]=0



      for (i=0;i<loc1.size();i++)
      {
            if (loc1[i]<loc2[0])
            {
                t2.insert(i);
                loc1toBin[i]=0;
            }
      }
      t1[4]=t2;
//    ls1.append(t1)
      ls1[0]=t1;
//    for i in range(len(loc2)-1):
    for (i=0;i<loc2.size()-1;i++)
    {
//        t1['loc2start']=loc2[i]
      ls1[i+1].loc2Start=loc2[i];
//        t1['loc2end']=loc2[i+1]
      ls1[i+1].loc2End=loc2[i+1];
//        t1['loc2Height']=val2[i]
      ls1[i+1].loc2Height=loc2[i+1];
    }
//        t1['otherPoints'],temp1=getOtherPoints(loc1,val1,loc2[i],loc2[i+1])
//        if len(t1['otherPoints'])>0:
//            t1['loc1Height']=max([x for x in val1 if x<t1['otherPoints'][0][1] ])
//        else:
//            t1['loc1Height']=getHeightAtPoint(loc1,val1,loc2[i])
//        for j in temp1:
//            loc1toBin[j]=len(ls1)
//        ls1.append(t1)
//    t1['loc2start']=loc2[-1]
      ls1[loc2.size()-1].loc2Start=loc2[-1];
//    t1['loc2end']=1000
      ls1[loc2.size()-1].loc2End=1000;
//    t1['loc2Height']=val2[-1]
      ls1[loc2.size()-1].loc2Height=val2[val2.size()-1];

// okay lets add the remaining fields to the data structure.
    int j=0;
    for (i=0;i<loc1.size();i++)
    {
        while (loc2[j]<loc1[i])
        {
            ls1[j].loc1IndexStart=i;
            j+=1;
            if (j==loc2.size())
            {break;}
            ls1[j].loc1IndexStart=i;
        }
        loc1toBin[i]=j;
        //this is annoying could do with a goto here!
        if (j==loc2.size())
        {break;}
    }
    ls1[ls1.size()-1].loc1IndexStart=i;
    ls1[ls1.size()-1].loc1IndexEnd=loc1.size()-1;
    for (k=i;k<loc1.size()-1;i++)
    {loc1toBin[k]=ls1.size()-1;}

//    ## first pass on the data.
//    res=0
      double res;
//    diffs=[]
    std::vector<double> diffs(ls1.size())
    int s1;
    int se;
    double v1;
    double h1l1;
    double prevH;
//    for item in ls1:
    for (i=0;i<ls1.size();i++)
    {
//        s1=item['loc2start']
        s1=ls1[i].loc2start;
//        se=item['loc2end']
        se=ls1[i].loc2end;
//        t1=item
//        v1=item['loc2Height']
        v1=ls1[i].loc2Height;
//        h1l1=item['loc1Height']
        h1l1=ls1[i].loc1Height;
//        if len(t1['otherPoints'])>0:
        if (ls1[i].loc1IndexEnd!=ls1[i].loc1IndexStart)
        {
//            prevH=h1l1
            prevH=h1l1;
//            for p1 in t1['otherPoints']:
            for (j=ls1[i].loc1IndexStart;j<ls1[i].loc1IndexEnd;j++)
            {
//                res+=(p1[0]-s1)*abs(v1-prevH)
                    res+=(loc1[j]-s1)*abs(v1-prevH);
//                prevH=p1[1]
                prevH=val1[j];
//                s1=p1[0]
                s1=loc1[j];
            }
//            # final point
//            res+=(se-f1[0])*abs(f1[1]-v1)
            res+=(se-loc1[ls1[i].loc1IndexEnd])*abs(val1[ls1[i].loc1IndexEnd]-v1) // could be an out by 1 error here.
//            diffs.append(-abs(f1[1]-v1)+abs(h1l1-v1))
            diffs[i]=abs(h1l1-v1)-abs(val1[ls1[i].loc1IndexEnd]-v1);
        }
        else
        {
//            res+=(item['loc2end']-item['loc2start'])*abs(item['loc2Height']-item['loc1Height'])
            res+=(ls1[i].loc2End-ls1[i].loc2Start)*abs(ls1[i].loc2Height-ls1[i].loc1Height);
//            diffs.append(0)
            diffs[i]=0;
        }
    }

    //completed first pass, now lets continue
    double curBestEmd=res;
    double curBestOffset=minOffset;
//    diffs1=sum(diffs)
     double diffs1=0;
     for (i=0;i<diffs.size();i++)
     {diffs1+=diffs[i];}
//    oldOffset=0
    double oldOffset=0
    double resOld;
int bin2Update;
//    while len(offsets)>0:
//        if item[0]==0:
//            continue
    for (i=0;i<offsets.size();i++)
    {
        if (offsets[i][0]!=0)
        {break;}
    }
    i-=1;
    while (i<offsets.size())
    {
        i+=1;
//        res+=(item[0]-oldOffset)*diffs1
//        r
        res+=(offsets[i].offset-oldOffset)*diffs1;
        if (res<curBestEmd)
        {
            curBestEmd=res;
            curBestOffset=offsets[i].offset;
        }
//        oldOffset=item[0]
        oldOffset=offsets[i].offset;
//        q1=[item,]
//        if len(offsets)>0:
//            while offsets[0][0]==item[0]:
//                q1.append(offsets.pop(0))
//        for item1 in q1:
        while (offsets[i].offset==oldOffset)
        {
//            bin2Update=loc1toBin[item1[1]]
            bin2Update=loc1toBin[offsets[i].loc1loc];
//            ## need to update this bin
//            loc1toBin[item1[1]]+=1
            loc1toBin[offsets[i].loc1loc]+=1;
//            aq1=ls1[bin2Update]
//            sd1=aq1['otherPoints'].pop()
            ls1[bin2Update].loc1IndexEnd-=1;
//            ls1[bin2Update+1]['otherPoints']=[sd1,]+ls1[bin2Update+1]['otherPoints']
            ls1[bin2Update+1].loc1IndexStart-=1;
//            t1=ls1[bin2Update+1]
//            if len(t1['otherPoints'])>0:
//                t1['loc1Height']=max([0,]+[x for x in val1 if x<t1['otherPoints'][0][1] ])
//            else:
//                t1['loc1Height']=getHeightAtPoint(loc1,val1,loc2[-1])
            ls1[bin2Update+1].loc1Height=val1[ls1[bin2Update+1].loc1IndexStart]; // i think that this right??
//            t1=ls1[bin2Update]
//            if len(t1['otherPoints'])>0:
            if (ls1[bin2Update].loc1IndexEnd!=ls1[bin2Update].loc1IndexStart)
            {
//                f1=t1['otherPoints'][-1]
//                f0=t1['otherPoints'][0]
//                h1l1=t1['loc1Height']
                  h1l1=ls1[bin2Update].loc1Height;
//                v1=t1['loc2Height']
                v1=ls1[bin2Update].loc2Height;
//                diff1=-abs(f1[1]-v1)+abs(h1l1-v1)
                diff1=abs(h1l1-v1)-abs(val1[ls1[bin2Update].loc1IndexEnd]-v1);
            }
//            else:
//                diff1=0
            else
            {diff1=0;}

//            diffs1+=diff1-diffs[bin2Update]
            diffs1+=diff1-diffs[bin2Update];
//            diffs[bin2Update]=diff1
            diffs[bin2Update]=diff1;

            i+=1;
        }
    }
    std::cout << bestEmd << " " <<bestOffset <<"\n";
    return bestEmd;
}
