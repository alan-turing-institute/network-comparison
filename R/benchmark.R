library(tictoc)
library(igraph)
source('cppInterfaces.R')

a1l=c()

a2l=c()
a3l=c()
t1l=c()
t2l=c()
t3l=c()

xs=(1:20)*100
for (i in xs)
{
    G1=erdos.renyi.game(i,0.05)
    dhist1=gdd(G1)
    G2=erdos.renyi.game(i,0.05)
    dhist2=gdd(G2)
    ts=Sys.time();
    a1=net_emd_fast(dhist1,dhist2,method='optimise');
    t1=Sys.time()-ts;
    ts=Sys.time();
    a2=net_emd_fast(dhist1,dhist2,method='exhaustive');
    t2=Sys.time()-ts;
    ts=Sys.time();
    a3=net_emd_fast(dhist1,dhist2,method='exhaustiveVer2');
    t3=Sys.time()-ts;
    a1l<-c(a1l,a1)
    a2l<-c(a2l,a2)
    a3l<-c(a3l,a3)
    t1l<-c(t1l,t1)
    t2l<-c(t2l,t2)
    t3l<-c(t3l,t3)
}
