#tab page 26
rweib=rweibull;rlnor=rlnorm;rgamm=rgamma;rlogi=rlogis;
pweib=pweibull;plnor=plnorm;pgamm=pgamma;plogi=plogis;
escalw=function(sig,resol=3,x){
  dd <- wd(sig, filter.number=1, family="DaubExPhase");
  esc=accessC(dd,resol)/sqrt(2)^(log2(length(sig))-resol);#length(esc)
  t=seq(min(x),max(x),(max(x)-min(x))/(2^(resol)) );
  t=t[1:2^resol];#length(t);(esc)
  iso=isoreg(t,esc); 
  stepfun(t,c(iso$y[1],iso$y));#rm(t)
}
KSfw1=function(refcdf=HDL[[6]],cdftheo=pnorm,p1=.01,p2=1){
  x1=knots(refcdf);x2=shift(x1,1);x2[1]=x1[1];
  DnP=abs(cdftheo(x1,p1,p2)-refcdf(x1));DnM=abs(cdftheo(x1,p1,p2)-refcdf(x2));
  max(DnP,DnM)
}
VMfw1=function(refcdf=HDL,cdftheo=pnorm,p1=.01,p2=1){
  x1=knots(refcdf);F0=cdftheo(x1,p1,p2);Fn=refcdf ;#hist(F0);#hist(Fn(x));
  F0=shift(F0,1);F0[1]=0;F0[length(F0)+1]=1;
  dF0=diff(F0);F0=cdftheo(x1,p1,p2);
  sqrt(sum( ((Fn(x1)-F0)^2)*dF0 ))
}
MEfw1=function(wvlet,CDF,p1=.01,p2=1){
  x=knots(wvlet)
  return( mean(abs(wvlet(x)-CDF(x,p1,p2))) )
}

shapeDist=function(pdist,parm1=1,parm2=1){
  x=seq(-200,200,.0001);#y=pweibull(x,5,1);plot(x,y,'l')
  y=pdist(x,parm1,parm2);
  maxseq<-2*(x[which(y>.98)[1]]);minseq<-2*(x[which(y>.0001)[1]]);rm(x)
  x<-seq(minseq,maxseq,(maxseq-minseq)/2^12);x<-x[1:(2^12)];
  #plot(x,pdist(x,parm1,parm2),'l')
  list(x,pdist(x,parm1,parm2));
}

HaarGenGr=function(cdffun=dist){ 
  x=cdffun[[1]];y=cdffun[[2]];
  HDL=lapply(seq(3,6,1),function(resol){escalw(y,resol,x)});HDL[[1]]
  plot(x,y,"l",col="red",xlab="",ylab="",lwd=2,yaxt='n');axis(2,c(.1,1));par(new=T);
  plot(lwd=2,x,HDL[[4]](x),ylim=c(0,1),ylab = "",yaxt='n',"l",xlab="");
}
HaarWGene=function(cdffun){ 
  x=cdffun[[1]];y=cdffun[[2]];
  HDL=lapply(seq(3,7,1),function(resol){escalw(y,resol,x)});
}# we put numbers here to reduce the compute time

# pas envie de boucle car je suis cerveau est boucle

# dist=shapeDist(plnor,1,1)    ; HaarGenGr(dist);
# dist=shapeDist(pnorm,0,1)    ; HaarGenGr(dist); #VMf(y[[4]],pnorm)
# dist=shapeDist(pweib,1,1)    ; HaarGenGr(dist);
# dist=shapeDist(plogi,2,1)    ; HaarGenGr(dist);
# dist=shapeDist(pgamm,1,2)    ; HaarGenGr(dist);
lev=5;#3==5 resolution;5 is 7 resolution
cat(c(KSfw1(HaarWGene(shapeDist(pnorm,0,1))[[lev]], pnorm,0,1),
      KSfw1(HaarWGene(shapeDist(plnor,1,1))[[lev]], plnor,1,1),
      KSfw1(HaarWGene(shapeDist(pweib,3,1))[[lev]], pweib,3,1),
      KSfw1(HaarWGene(shapeDist(pgamm,1,2))[[lev]], pgamm,1,2),
      KSfw1(HaarWGene(shapeDist(plogi,2,1))[[lev]], plogi,2,1)),"\n");

cat(c(VMfw1(HaarWGene(shapeDist(pnorm,0,1))[[lev]], pnorm,0,1),
      VMfw1(HaarWGene(shapeDist(plnor,1,1))[[lev]], plnor,1,1),
      VMfw1(HaarWGene(shapeDist(pweib,3,1))[[lev]], pweib,3,1),
      VMfw1(HaarWGene(shapeDist(pgamm,1,2))[[lev]], pgamm,1,2),
      VMfw1(HaarWGene(shapeDist(plogi,2,1))[[lev]], plogi,2,1)),"\n");

cat(c(MEfw1(HaarWGene(shapeDist(pnorm,0,1))[[lev]], pnorm,0,1),
      MEfw1(HaarWGene(shapeDist(plnor,1,1))[[lev]], plnor,1,1),
      MEfw1(HaarWGene(shapeDist(pweib,3,1))[[lev]], pweib,3,1),
      MEfw1(HaarWGene(shapeDist(pgamm,1,2))[[lev]], pgamm,1,2),
      MEfw1(HaarWGene(shapeDist(plogi,2,1))[[lev]], plogi,2,1)),"\n")

vmfnorm=NULL;vmflnor=NULL;vmfweib=NULL;vmfgamm=NULL;vmflogi=NULL;
KSfnorm=NULL;KSflnor=NULL;KSfweib=NULL;KSfgamm=NULL;KSflogi=NULL


###########################################################3
tmp=NULL;lev=2^7;x=NULL;
vmfnorm=NULL;vmflnor=NULL;vmfweib=NULL;vmfgamm=NULL;vmflogi=NULL;
KSfnorm=NULL;KSflnor=NULL;KSfweib=NULL;KSfgamm=NULL;KSflogi=NULL;
MEfnorm=NULL;MEflnor=NULL;MEfweib=NULL;MEfgamm=NULL;MEflogi=NULL

VMfs1=function(sampdist,cdftheo=pnorm,p1=0,p2=1){
  if(!is.numeric(x))x<<-shapeDist(cdftheo,p1,p2)[[1]];
  F0=cdftheo(x,p1,p2);Fn=sampdist;
  F0=shift(F0,1);F0[1]=0;F0[length(F0)+1]=1;# cdf limits.
  dF0=diff(F0);F0=cdftheo(x,p1,p2);
  sqrt(sum( ((Fn(x)-F0)^2)*dF0 ))
}
KSfs1=function(sampdist=ecdf(rnorm(lev,0,1)),cdftheo=pnorm,p1=0,p2=1){
  if(!exists("x1"))x1<<-shapeDist(cdftheo,p1,p2)[[1]];
  x2=shift(x1,1);x2[1]=x1[1];
  DnP=abs(cdftheo(x1,p1,p2)-sampdist(x1));
  DnM=abs(cdftheo(x1,p1,p2)-sampdist(x2));
  max(DnP,DnM)
}
MEfs1=function(sampdist=ecdf(rnorm(lev,0,1)),cdftheo=pnorm,p1=0,p2=1){
  if(!exists("x2"))x2=shapeDist(cdftheo,p1,p2)[[1]];
  return( mean(abs(sampdist(x2)-cdftheo(x2,p1,p2)) ) ) 
}
p=1024;
for(i in 1:p) KSfnorm=(c(KSfnorm,KSfs1(ecdf(rnorm(lev,0,1)),pnorm,0,1)));
for(i in 1:p) KSflnor=(c(KSflnor,KSfs1(ecdf(rlnor(lev,1,1)),plnor,1,1)));
for(i in 1:p) KSfweib=(c(KSfweib,KSfs1(ecdf(rweib(lev,3,1)),pweib,3,1)));
for(i in 1:p) KSfgamm=(c(KSfgamm,KSfs1(ecdf(rgamm(lev,1,2)),pweib,1,2)));
for(i in 1:p) KSflogi=(c(KSflogi,KSfs1(ecdf(rlogi(lev,2,1)),plogi,2,1)))
cat(c(median(KSfnorm),median(KSflnor),median(KSfweib),median(KSfgamm),median(KSflogi)),"\n");
cat( c(sum(KSfnorm>.01),sum(KSflnor>.04),sum(KSfweib>.01),sum(KSfgamm >.03),sum(KSflogi>.02))/p,'\n');

for(i in 1:p) vmfnorm=(c(vmfnorm,VMfs1(ecdf(rnorm(lev,0,1)),pnorm,0,1)));
for(i in 1:p) vmflnor=(c(vmflnor,VMfs1(ecdf(rlnor(lev,1,1)),plnor,1,1)));
for(i in 1:p) vmfweib=(c(vmfweib,VMfs1(ecdf(rweib(lev,3,1)),pweib,3,1)));
for(i in 1:p) vmfgamm=(c(vmfgamm,VMfs1(ecdf(rgamm(lev,1,2)),pweib,1,2)));
for(i in 1:p) vmflogi=(c(vmflogi,VMfs1(ecdf(rlogi(lev,2,1)),plogi,2,1)));
cat(c(median(vmfnorm),median(vmflnor),median(vmfweib),median(vmfgamm),median(vmflogi)),"\n");
cat( c(sum(vmfnorm>.01),sum(vmflnor>.02),sum(vmfweib>.01),sum(vmfgamm >.01),sum(vmflogi>.01))/p,'\n' );

for(i in 1:p) MEfnorm=(c(MEfnorm,MEfs1(ecdf(rnorm(lev,0,1)),pnorm,0,1)));
for(i in 1:p) MEflnor=(c(MEflnor,MEfs1(ecdf(rlnor(lev,1,1)),plnor,1,1)));
for(i in 1:p) MEfweib=(c(MEfweib,MEfs1(ecdf(rweib(lev,3,1)),pweib,3,1)));
for(i in 1:p) MEfgamm=(c(MEfgamm,MEfs1(ecdf(rgamm(lev,1,2)),pweib,1,2)));
for(i in 1:p) MEflogi=(c(MEflogi,MEfs1(ecdf(rlogi(lev,2,1)),plogi,2,1)))
cat(c(median(MEfnorm),median(MEflnor),median(MEfweib),median(MEfgamm),median(MEflogi)),"\n")
vmfnorm=NULL;vmflnor=NULL;vmfweib=NULL;vmfgamm=NULL;vmflogi=NULL;
KSfnorm=NULL;KSflnor=NULL;KSfweib=NULL;KSfgamm=NULL;KSflogi=NULL



