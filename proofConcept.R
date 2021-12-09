n=2^11;xseq=seq(-3.9999,4, 8/n);ycdf=pnorm(xseq);ns=128;
graph=function(rmin=2,rmax=11){
  par(mar=c(5,5,4,4)+.4);
  plot(rmin:rmax,(wNdistKS[-c( 1:(rmin-1),(rmax+1):(log2(n)+1))]),"l",
       lwd=3,xlab="",ylab="",col="blue",font=1,col.axis="blue",
       xlim=c(rmin,rmax),ylim=c(0,80),cex.axis=1.4);
  par(new=T);
  axis(2, font=2,col.axis = "blue",cex.axis=1.4);
  par(new=T);
  plot(rmin:rmax,(Byte[-c( 1:(rmin-1),(rmax+1):(log2(n)+1))]),lwd=3,
       ylab = "",xlab = "",col="red","l",
       font=2,cex.main=1.,ylim=c(0,80),yaxt='n',xaxt='n');
  mtext(side=2,line=3,expression(bold('%  KS Error ')),font=6,col="blue",ylim=c(0,50),cex=1.5);
  mtext(side=4,line=3,expression(bold('%  Byte rate')),font=6,col="red",cex=1.5);
  mtext(side=1,line=3,expression(bold(Resolution)),font=6,cex=1.6)
  mtext(side=3,line=1,expression(bold("Save Memory, limit Error")),font=6,cex=1.8)
  axis(4,at=seq(0,80,10),font=2,col.axis="red",cex.axis=1.4)
  axis(1,at=seq(0,12,2),font=2,col.axis="black",cex.axis=1.4)
  abline(v=7.2, col="chartreuse", lwd=3, lty=2)
  abline(h=3.3, col="chartreuse",lwd=3, lty=2)
  rect(6, 0, 8, 10,border = "orange2",lwd = 3)
}
## escal refer to escalier. in fact Haar is stairs
escal=function(sig,resol=3){
  dd <- wd(sig, filter.number=1, family="DaubExPhase");
  esc=accessC(dd,resol)/sqrt(2)^(log2(n)-resol);#length(esc)
  t=seq( -3.9999999,4,8/(2^(resol)) );#length(t);length(esc)
  iso=isoreg(t,esc); 
  stepfun(t,c(iso$y[1],iso$y));
}

####extractSize of object from wavelet############
extractSizW=function(HDL){
  return(as.integer(object.size(knots(HDL)))*100/as.integer(object.size(ycdf)))
}
Normalise=function(v){
  v=abs(v-v[log2(n)+1])*100/abs(v[1]-v[log2(n)+1]);
  ifelse(v>100,100,v);
}
###stat distances 
KSf=function(refcdf=HDL[[6]],cdftheo=pnorm){
  x1=knots(refcdf);x2=shift(x1,1);x2[1]=x1[1];
  DnP=abs(cdftheo(x1)-refcdf(x1));DnM=abs(cdftheo(x1)-refcdf(x2));
  max(DnP,DnM)
};


HDL=lapply(seq(0,log2(n),1),function(resol){escal(ycdf,resol)});
###############################################################
Byte=unlist(lapply(HDL, extractSizW));
######distances to wavlet #####################################
wdistKS =unlist(lapply(HDL,KSf,pnorm));  wNdistKS =Normalise(wdistKS);


graph()
