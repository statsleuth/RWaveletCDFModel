#set.seed(364575). pnorm committed error while approacing with 7 resolution wavelet
pck=c("goftest","data.table","latex2exp","gridExtra","MASS","mvtnorm",
      "Morpho","wavethresh","Metrics","crayon")
#install.packages(pck);
lapply(pck[], library, character.only = TRUE,quietly=TRUE)
mr=2^12;lev=7;#maxreso;
escalw=function(sig,resol=3,x){
  dd <- wd(sig, filter.number=1, family= "DaubExPhase" );
  esc=accessC(dd,resol)/sqrt(2)^(log2(length(sig))-resol);#length(esc)
  t=seq(min(x),max(x),(max(x)-min(x))/(2^(resol)) );
  t=t[1:2^resol];#length(t);(esc)
  iso=isoreg(t,esc); 
  stepfun(t,c(iso$y[1],iso$y));#rm(t)
}
shapeDist=function(pdist,parm1=1,parm2=1){
  x=seq(-200,200,.0001);y=pdist(x,parm1,parm2);
  maxseq<-2*(x[which(y>.98)[1]]);minseq<-2*(x[which(y>.001)[1]]);
  x=NULL;x<-seq(minseq,maxseq,(maxseq-minseq)/mr);x<-x[1:mr];
  list(x,pdist(x,parm1,parm2));
}
HaarWGene=function(cdffun){ 
  x=cdffun[[1]];y=cdffun[[2]];
  HDL=lapply(seq(1,7,1),function(resol){escalw(y,resol,x)});
}
KSfw=function(refcdf=HDL[[6]],cdftheo=pnorm,p1=.01,p2=1){
  x1=knots(refcdf);x2=shift(x1,1);x2[1]=x1[1];
  DnP=abs(cdftheo(x1,p1,p2)-refcdf(x1));DnM=abs(cdftheo(x1,p1,p2)-refcdf(x2));
  max(DnP,DnM)
}
VMfw=function(refcdf=HDL,cdftheo=pnorm,p1=.01,p2=1){
  x1=knots(refcdf);F0=cdftheo(x1,p1,p2);Fn=refcdf ;#hist(F0);#hist(Fn(x));
  F0=shift(F0,1);F0[1]=0;F0[length(F0)+1]=1;
  dF0=diff(F0);F0=cdftheo(x1,p1,p2);
  sqrt(sum( ((Fn(x1)-F0)^2)*dF0 ))
}
lev=7
HDL=HaarWGene(shapeDist(pnorm,0,1))[[lev]];
KSfw(HDL,pnorm,0,1);
VMfw(HDL, pnorm,0,1)
