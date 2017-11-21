

def getHeightAtPoint(loc1,val1,point):
    res=0
    for i in range(len(loc1)):
        if point>=loc1[i]:
            res=val1[i]
        else:
            return res
    return res

def getOtherPoints(loc1,val1,start1,end1):
    res=[]
    res1=[]
    for i in range(len(loc1)):
        if start1<=loc1[i]:
            if end1>loc1[i]:
                res.append([loc1[i],val1[i]])
                res1.append(i)
    print res1
    return res,res1

def constantVersionWithAccumPythonVersion(loc1,val1,loc2,val2,offsets):
    offsets=[]
    for i in range(len(loc1)):
        for j in range(len(loc2)):
            offsets.append([loc2[j]-loc1[i],i,j])
    offsets.sort()
    for i in range(len(loc1)):
        loc1[i]+=offsets[0][0]
    temp1=offsets[0][0]
    minOffset=offsets[0][0]
    for item in offsets:
        item[0]-=temp1
    loc1toBin={}
    ls1=[]
    t1={}
    t1['loc2start']=-10000
    t1['loc2end']=loc2[0]
    t1['loc2Height']=0
    t1['loc1Height']=getHeightAtPoint(loc1,val1,-10000)
    t1['otherPoints'],temp1=getOtherPoints(loc1,val1,-10000,loc2[0])
    for j in temp1:
        loc1toBin[j]=0
    ls1.append(t1)
    for i in range(len(loc2)-1):
        t1={}
        t1['loc2start']=loc2[i]
        t1['loc2end']=loc2[i+1]
        t1['loc2Height']=val2[i]
        t1['otherPoints'],temp1=getOtherPoints(loc1,val1,loc2[i],loc2[i+1])
        if len(t1['otherPoints'])>0:
            t1['loc1Height']=max([x for x in val1 if x<t1['otherPoints'][0][1] ])
        else:
            t1['loc1Height']=getHeightAtPoint(loc1,val1,loc2[i])
        for j in temp1:
            loc1toBin[j]=len(ls1)
        ls1.append(t1)
    t1={}
    t1['loc2start']=loc2[-1]
    t1['loc2end']=1000
    t1['loc2Height']=val2[-1]
    t1['otherPoints'],temp1=getOtherPoints(loc1,val1,loc2[-1],10000)
    if len(t1['otherPoints'])>0:
        t1['loc1Height']=max([x for x in val1 if x<t1['otherPoints'][0][1] ])
    else:
        t1['loc1Height']=getHeightAtPoint(loc1,val1,loc2[-1])
    ls1.append(t1)
    for j in temp1:
        loc1toBin[j]=len(ls1)
    ls1.append(t1)
    ## first pass on the data.
    res=0
    diffs=[]
    count=0
    for item in ls1:
        print 'res',count,res
        count+=1
        t1=item
        s1=item['loc2start']
        se=item['loc2end']
        v1=item['loc2Height']
        h1l1=item['loc1Height']
        print '2',count,res
        if len(t1['otherPoints'])>0:
            prevH=h1l1
            for p1 in t1['otherPoints']:
                res+=(p1[0]-s1)*abs(v1-prevH)
                print 'res update mark2 ',res
                prevH=p1[1]
                s1=p1[0]
            f1=t1['otherPoints'][-1]
            f0=t1['otherPoints'][0]
            # final point
            res+=(se-f1[0])*abs(f1[1]-v1)
            print '2.5',count,res
            diffs.append(-abs(f1[1]-v1)+abs(h1l1-v1))
        else:
            res+=(item['loc2end']-item['loc2start'])*abs(item['loc2Height']-item['loc1Height'])
            print '3.5',count,res
            ## this is wrong
            diffs.append(0)
        print '3',res
    results=[]
    results.append([res,minOffset])
    diffs1=sum(diffs)
    print diffs
    oldOffset=0
    while len(offsets)>0:
        item=offsets.pop(0)
        if item[0]==0:
            continue
        print 'diffs=',diffs
        ## okay lets take the first offset
        res+=(item[0]-oldOffset)*diffs1
        results.append([res,minOffset+item[0],diffs1])
        oldOffset=item[0]
        q1=[item,]
        if len(offsets)>0:
            while offsets[0][0]==item[0]:
                q1.append(offsets.pop(0))
        for item1 in q1:
            ## need to update this bin
            bin2Update=loc1toBin[item1[1]]
            loc1toBin[item1[1]]+=1
            aq1=ls1[bin2Update]
            sd1=aq1['otherPoints'].pop()
            ls1[bin2Update+1]['otherPoints']=[sd1,]+ls1[bin2Update+1]['otherPoints']
            t1=ls1[bin2Update+1]
            if len(t1['otherPoints'])>0:
                t1['loc1Height']=max([0,]+[x for x in val1 if x<t1['otherPoints'][0][1] ])
            else:
                t1['loc1Height']=getHeightAtPoint(loc1,val1,loc2[-1])
            t1=aq1
            t1=ls1[bin2Update]
            if len(t1['otherPoints'])>0:
                f1=t1['otherPoints'][-1]
                f0=t1['otherPoints'][0]
                h1l1=t1['loc1Height']
                v1=t1['loc2Height']
                diff1=-abs(f1[1]-v1)+abs(h1l1-v1)
                #diff1=abs(f0[1]-v1)-abs(f1[1]-v1)
            else:
                diff1=0
            diffs1+=diff1-diffs[bin2Update]
            diffs[bin2Update]=diff1
            bin2Update+=1
            t1=ls1[bin2Update]
            if len(t1['otherPoints'])>0:
                f1=t1['otherPoints'][-1]
                f0=t1['otherPoints'][0]
                v1=t1['loc2Height']
                h1l1=t1['loc1Height']
                diff1=-abs(f1[1]-v1)+abs(h1l1-v1)
            else:
                diff1=0
            diffs1+=diff1-diffs[bin2Update]
            diffs[bin2Update]=diff1
    return results


def constantVersionWithAccum(loc1,val1,loc2,val2):
    currentEmd=None
    loc1=[min(loc1[0],loc2[0])-1,]+loc1
    loc2=[min(loc1[0],loc2[0])-1,]+loc2
    val1=[0,]+val1
    val2=[0,]+val2
    offsets=[]
    for i in range(len(loc1)):
        for j in range(len(loc2)):
            offsets.append([loc2[j]-loc1[i],i,j])
    offsets.sort()

    ## updating offsets to make everything positive (dont have to do this but
    ##makes life easier to think about)
    for i in range(len(loc1)):
        loc1[i]+=offsets[0][0]
    temp1=offsets[0][0]
    for item in offsets:
        item[0]+=temp1

    interval=[]
    vals=[]
    ## first pass
    for i in range(1,len(loc1)):
        notFound=True
        for j in range(1,len(loc2)):
            if loc1[i]<loc2[j]:
                interval.append(j)
                h1=abs(val1[i-1]-val2[j-1])
                h2=abs(val1[i]-val2[j-1])
                width=loc1[i]-loc2[j-1]
                vals.append(h2-h1)
                notFound=False
                break
        if notFound:
            inverval.append(len(loc2))
    sum(vals)
