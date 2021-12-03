# execute apprCDF1 first

vmfnorm=NULL;vmflnor=NULL;vmfweib=NULL;vmfgamm=NULL;vmflogi=NULL;
KSfnorm=NULL;KSflnor=NULL;KSfweib=NULL;KSfgamm=NULL;KSflogi=NULL;

escalw=function(sig,resol=3,x){
  dd <- wd(sig, filter.number=1, family= "DaubExPhase" );
  esc=accessC(dd,resol)/sqrt(2)^(log2(length(sig))-resol);#length(esc)
  t=seq(min(x),max(x),(max(x)-min(x))/(2^(resol)) );
  t=t[1:2^resol];#length(t);(esc)
  iso=isoreg(t,esc); 
  stepfun(t,c(iso$y[1],iso$y));#rm(t)
}
shapeComp=function(dist1,dist2,p1=3,p2=4.2,p3=4.5,p4=3.2){
  x1=shapeDist(dist1,p1,p2)[[1]];x2=shapeDist(dist2,p3,p4)[[1]];
  x1=round(max(abs(min(x1,x2)),abs(max(x1,x2))))
  xseq=seq(-x1,x1,2*x1/2^10);xseq=xseq[1:2^10]
  list(xseq,.4*dist1(xseq,p1,p2)+.6*dist2(xseq,p3,p4));
}
#shapeComp(pweibull,pweibull,3,4.2,4.5,3.2);
KSfw=function(refcdf=HDL[[6]],cdfcomp){
  x1=cdfcomp[[1]];x2=shift(x1,1);x2[1]=x1[1];
  DnP=abs(cdfcomp[[2]]-refcdf(x1));DnM=abs(cdfcomp[[2]]-refcdf(x2));
  max(DnP,DnM)
}
VMfw=function(refcdf=HaarWGene(c1)[[4]],cdfcomp=c1){
  x1=cdfcomp[[1]];F0=cdfcomp[[2]];Fn=refcdf ;#hist(F0);#hist(Fn(x));
  F0=shift(F0,1);F0[1]=0;F0[length(F0)+1]=1;
  dF0=diff(F0);F0=cdfcomp[[2]];
  sqrt(sum( ((Fn(x1)-F0)^2)*dF0 ))
}
MEfw=function(wvlet=HaarWGene(mix)[[1]],cdfcomp=mix){
  x=cdfcomp[[1]]
  return( mean(abs(wvlet(x)-cdfcomp[[2]])) )
}

mix=shapeComp(plnor,plnor,1,1,3.3,.08);  cat(KSfw(HaarWGene(mix)[[6]],mix),"\n");HaarGenGr(mix);
mix=shapeComp(pweib,pweib,3,4.2,4.5,5);  cat(KSfw(HaarWGene(mix)[[6]],mix),"\n")

mix=shapeComp(plnor,plnor,1,1,3.3,.08);  cat(VMfw(HaarWGene(mix)[[6]],mix),"\n")
mix=shapeComp(pweib,pweib,3,4.2,4.5,5);  cat(VMfw(HaarWGene(mix)[[6]],mix),"\n")

mix=shapeComp(plnor,plnor,1,1,3.3,.08);  cat(MEfw(HaarWGene(mix)[[6]],mix),"\n")
mix=shapeComp(pweib,pweib,3,4.2,4.5,5);  cat(MEfw(HaarWGene(mix)[[6]],mix),"\n")





