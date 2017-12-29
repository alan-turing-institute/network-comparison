library(igraph)
source('emd.R')
source('measures_net_emd.R')
source('dhist.R')
source('orca_interface.R')
source('cppInterfaces.R')

a0l=c()
a1l=c()
a2l=c()
a3l=c()
t0l=c()
t1l=c()
t2l=c()
t3l=c()
xs=(1:35)*100
xs=(1:30)*500

for (i in xs)
{
    G1=erdos.renyi.game(i,0.05)
    dhist1=gdd(G1)
    G2=barabasi.game(i,5)
    dhist2=gdd(G2)
    ts=Sys.time();
    a0=net_emd(dhist1,dhist2,method='optimise');
    t0=Sys.time()-ts;
    ts=Sys.time();
    a1=net_emd_fast(dhist1,dhist2,method='optimise');
    t1=Sys.time()-ts;
    ts=Sys.time();
    a2=net_emd_fast(dhist1,dhist2,method='exhaustive');
    t2=Sys.time()-ts;
    ts=Sys.time();
#    a3=net_emd_fast(dhist1,dhist2,method='exhaustiveVer2');
#    t3=Sys.time()-ts;
    a3=0
    t3=0
    a0l<-c(a0l,a0)
    a1l<-c(a1l,a1)
    a2l<-c(a2l,a2)
    a3l<-c(a3l,a3)
    t0l<-c(t0l,t0)
    t1l<-c(t1l,t1)
    t2l<-c(t2l,t2)
    t3l<-c(t3l,t3)
}


plot(xs,t2l,log="y",col='blue',ylim=c(0.0001,100),lty = 1)
lines(xs,t2l,col='blue')
points(xs,t1l,col='red',lty = 2)
lines(xs,t1l,col='red')
points(xs,t3l,col='green',lty = 3)
lines(xs,t3l,col='green')
points(xs,t0l,col='black',lty = 4)
lines(xs,t0l,col='black')
legend("topleft",c("ExhaustfastVer1 ","OptFast","ExhaustFastVer2","Opt"),pch=c(1,1,1,1),  col=c("blue","red","green","black"))


