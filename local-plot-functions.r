#--plots for srac2lrac analysis
#--this R code by Rohan Williams (rbhwilliams@gmail.com)

plot.segments<-function(segmentData,...)
{
 segments(x0=segmentData$x0,y0=segmentData$y0,x1=segmentData$x1,y1=segmentData$y1,...)
}

summary.plot.for.srac2lrac.single<-function(srac2lracTableList,lracDNAStringSet,lrac.id,srac.id)
{
 srac2lracTableS<-srac2lracTableList[[lrac.id]][[srac.id]]
 if(nrow(srac2lracTableS)>0)
 {
  lrac.fasta<-lracDNAStringSet[lrac.id]
  #--turn this into "segment data"
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=rep(1,nrow(srac2lracTableS)),x1=srac2lracTableS[,"send"],y1=rep(1,nrow(srac2lracTableS)))
  xad<-rowMeans(cbind(srac2lracTableS[,"send"],srac2lracTableS[,"sstart"]))
 
  lo<-layout(matrix(c(1,2,3,4,5),5,1,byrow=T),height=c(0.8,rep(0.7,3),1),width=rep(1,5))
  par(mar=c(0.1,6.1,1.5,2.1))
  plot(xad,rep(1,nrow(srac2lracTableS)),type="n",xaxt="n",yaxt="n",las=1,ylab="",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  legend(0,1.4,srac.id,ncol=1,pch=16,col=1,cex=1,bty="n")
  mtext(side=2,at=1.4,line=4,"A.",cex=1,las=1)
  par(mar=c(0.1,6.1,0.1,2.1))
  plot(xad,srac2lracTableS[,"pident"],xaxt="n",las=1,ylab="PID",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
  mtext(side=2,at=100,line=4,"B.",cex=1,las=1)
  par(mar=c(0.1,6.1,0.1,2.1))
  plot(xad,srac2lracTableS[,"al2ql"],xaxt="n",xlab="",las=1,ylab="al2ql",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1.12))
  mtext(side=2,at=1.1,line=4,"C.",cex=1,las=1)
  par(mar=c(0.1,6.1,0.1,2.1))
  plot(xad,log10(srac2lracTableS[,"bitscore"]),xaxt="n",xlab="",las=1,ylab="Bitscore (log10)",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(1.5,6.5))
  mtext(side=2,at=6.5,line=4,"D.",cex=1,las=1)
  par(mar=c(5.1,6.1,0.1,2.1))
  plot(xad,srac2lracTableS[,"qgc"],las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="Query GC",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1))
  mtext(side=2,at=6.5,line=4,"E.",cex=1,las=1)
  par(mfrow=c(1,1))
 }
 else
 if(nrow(srac2lracTableS)<1)
 {
  print("No alignments are available for this pair!\n")
 }
 data.frame(srac2lracTableS,xad=xad,stringsAsFactors=F)
}

summary.plot.for.srac2lrac.single.aug<-function(srac2lracAugRes,lracDNAStringSet)
{
 #--'srac2lracResAug' is an annotated augmented element of an 'srac2lracTableList'
 
 lrac.fasta<-lracDNAStringSet
 #--turn this into "segment data" (we will reuse this variable multiple times
 srac2lracAugRes.segmentData<-data.frame(x0=srac2lracAugRes[,"sstart"],y0=rep(1,nrow(srac2lracAugRes)),x1=srac2lracAugRes[,"send"],y1=rep(1,nrow(srac2lracAugRes)))
 xad<-rowMeans(cbind(srac2lracAugRes[,"send"],srac2lracAugRes[,"sstart"]))
 
 lo<-layout(matrix(c(1,2,3,4,5),5,1,byrow=T),height=c(0.8,rep(0.7,3),1),width=rep(1,5))
 par(mar=c(0.1,6.1,1.5,2.1))
 plot(xad,rep(1,nrow(srac2lracAugRes)),type="n",xaxt="n",yaxt="n",las=1,ylab="",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)))
 plot.segments(srac2lracAugRes.segmentData,col=1,lwd=2)
 mtext(side=2,at=1.4,line=4,"A.",cex=1,las=1)
 par(mar=c(0.1,6.1,0.1,2.1))
 plot(xad,srac2lracAugRes[,"pident"],type="n",xaxt="n",las=1,ylab="PID",pch=3,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
 abline(v=xad,lty=2,col="lightgrey")
 srac2lracAugRes.segmentData<-data.frame(x0=srac2lracAugRes[,"sstart"],y0=srac2lracAugRes[,"pident"],x1=srac2lracAugRes[,"send"],y1=srac2lracAugRes[,"pident"])
 plot.segments(srac2lracAugRes.segmentData,col=1,lwd=2)
 mtext(side=2,at=100,line=4,"B.",cex=1,las=1)
 par(mar=c(0.1,6.1,0.1,2.1))
 plot(xad,srac2lracAugRes[,"al2ql"],type="n",xaxt="n",xlab="",las=1,ylab="al2ql",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1.12))
 abline(v=xad,lty=2,col="lightgrey")
 srac2lracAugRes.segmentData<-data.frame(x0=srac2lracAugRes[,"sstart"],y0=srac2lracAugRes[,"al2ql"],x1=srac2lracAugRes[,"send"],y1=srac2lracAugRes[,"al2ql"])
 plot.segments(srac2lracAugRes.segmentData,col=1,lwd=2)
 mtext(side=2,at=1.1,line=4,"C.",cex=1,las=1)
 par(mar=c(0.1,6.1,0.1,2.1))
 plot(xad,srac2lracAugRes[,"J.toi"],type="n",xaxt="n",xlab="",las=1,ylab="J-statistic",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
 abline(v=xad,lty=2,col="lightgrey")
 points(xad,srac2lracAugRes[,"J.toi"],pch=16,cex=0.90,col=1)
 points(xad,srac2lracAugRes[,"J.omx"],pch=21,bg="white",cex=0.90,col=1)
 mtext(side=2,at=100,line=4,"D.",cex=1,las=1)
 par(mar=c(5.1,6.1,0.1,2.1))
 plot(xad,log10(srac2lracAugRes[,"length"]),type="n",las=1,xlab=paste("Physical location on LR-chr (bp)",sep=" "),ylab="Query length (log10 bp)",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(1.3,6.3))
 abline(v=xad,lty=2,col="lightgrey")
 points(xad,log10(srac2lracAugRes[,"length"]),pch=3,cex=0.90,col=1)
 mtext(side=2,at=1.0,line=4,"E.",cex=1,las=1)
 par(mfrow=c(1,1))

 return(data.frame(srac2lracAugRes,xad=xad,stringsAsFactors=F))
}

#--added on 29 October 2019

#> par()$mar
#[1] 5.1 6.1 0.1 2.1

#--add some panels to put alignment data in context
summary.plot.for.srac2lrac.single<-function(srac2lracTableList,lracDNAStringSet,lrac.id,srac.id,blastn6.aug.df,contigSummData,sracBinMembershipList,Xlim,Ylim)
{
 srac2lracTableS<-srac2lracTableList[[lrac.id]][[srac.id]]
 if(nrow(srac2lracTableS)>0)
 {
  lrac.fasta<-lracDNAStringSet[lrac.id]
  #--turn this into "segment data"
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=rep(1,nrow(srac2lracTableS)),x1=srac2lracTableS[,"send"],y1=rep(1,nrow(srac2lracTableS)))
  xad<-rowMeans(cbind(srac2lracTableS[,"send"],srac2lracTableS[,"sstart"]))
 
  lo<-layout(matrix(c(1,2,3,4,5,6,6,6,7,7),5,2,byrow=F),height=c(0.8,rep(0.7,3),1),width=rep(1,5))
  par(mar=c(0.1,6.1,1.5,2.1))
  plot(xad,rep(1,nrow(srac2lracTableS)),type="n",xaxt="n",yaxt="n",las=1,ylab="",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  legend(0,1.4,srac.id,ncol=1,pch=16,col=1,cex=1,bty="n")
  mtext(side=2,at=1.4,line=4,"A.",cex=1,las=1)
  par(mar=c(0.1,6.1,0.1,2.1))
  plot(xad,srac2lracTableS[,"pident"],xaxt="n",las=1,ylab="PID",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
  mtext(side=2,at=100,line=4,"B.",cex=1,las=1)
  par(mar=c(0.1,6.1,0.1,2.1))
  plot(xad,srac2lracTableS[,"al2ql"],xaxt="n",xlab="",las=1,ylab="al2ql",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1.12))
  mtext(side=2,at=1.1,line=4,"C.",cex=1,las=1)
  par(mar=c(0.1,6.1,0.1,2.1))
  plot(xad,log10(srac2lracTableS[,"bitscore"]),xaxt="n",xlab="",las=1,ylab="Bitscore (log10)",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(1.5,6.5))
  mtext(side=2,at=6.5,line=4,"D.",cex=1,las=1)
  par(mar=c(5.1,6.1,0.1,2.1))
  plot(xad,srac2lracTableS[,"qgc"],las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="Query GC",pch=16,cex=0.90,col=1,cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1))
  mtext(side=2,at=1,line=4,"E.",cex=1,las=1)
  
  contigSummData.hb<-hexbin(contigSummData[,c("lgcov","gc")],xbins=100)
  
  par(mar=c(5.1,6.1,1.5,2.1))
  plot(contigSummData[,c("lgcov","gc")],pch=16,cex=0.75,col="lightgrey",las=1,xlab="Log coverage",ylab="GC",xlim=Xlim,ylim=Ylim)
  for(curbin in names(sracBinMembershipList))
  {
   Plot_ConvexHull(contigSummData[sracBinMembershipList[[curbin]],"lgcov"],contigSummData[sracBinMembershipList[[curbin]],"gc"],alpha("lightgrey",0.25),3)
  }
  points(contigSummData[sracBinMembershipList[[srac.id]],c("lgcov","gc")],pch=16,cex=0.75,col=1)
  legend(Xlim[1],Ylim[2],srac.id,ncol=1,pch=16,col=1,cex=1,bty="n")
  mtext(side=2,at=(Ylim[2]*1.01),line=4,"F.",cex=1,las=1)
  
  binNum<-substring(srac.id,5)
  print(binNum)
  kappaVal<-blastn6.aug.df$kappa[intersect(which(blastn6.aug.df$lrac==lrac.id),which(blastn6.aug.df$sracBin==binNum))]
  hd<-hist(blastn6.aug.df$kappa[which(blastn6.aug.df$lrac==lrac.id)],100,plot=F)
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count)*1.15),las=1,xlab=paste("Kappa statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1))
  arrows(kappaVal,max(hd$count)*1.10,kappaVal,max(hd$count),length=0.1)
  text(kappaVal,max(hd$count)*1.15,srac.id)
  mtext(side=2,at=max(hd$count)*1.15,line=4,"G.",cex=1,las=1)
  
 }
 else
 if(nrow(srac2lracTableS)<1)
 {
  print("No alignments are available for this pair!\n")
 }
 par(mfrow=c(1,1))
 data.frame(srac2lracTableS,xad=xad,stringsAsFactors=F)
}

Plot_ConvexHull<-function(xcoord,ycoord,lcolor,llty)
{
 hpts<-chull(x=xcoord,y=ycoord)
 hpts<- c(hpts,hpts[1])
 polygon(xcoord[hpts],ycoord[hpts],col=lcolor,lty=llty)
}

#--we need a function to compute GC content across an LR
calc.gc.content.across.lrac<-function(lrac,d=0.001)
{
 q<-round(quantile(1:nchar(lrac),probs=seq(0,1,d)),0)
 bw<-diff(q)
 cat("Bin widths will be between",min(bw),"and",max(bw),"bp in length\n")
 
 res=NULL
 gc=NULL
 mid=NULL
 for(i in (1:(length(q)-1)))
 {
  curmid<-mean(q[(i+1)],q[i])
  curseq<-subseq(lrac,q[i],q[(i+1)])
  curres<-oligonucleotideFrequency(curseq,1)[1,]
  curgc<-sum(curres[c("G","C")])/sum(curres)
  gc<-c(gc,curgc)
  mid<-c(mid,curmid)
 }
 res<-data.frame(mid=mid,gc=gc)
 res
}

set.zero.values.to.NA.in.hist.counts<-function(hdd)
{
 #--hdd is an output from 'hist'
 #--set values of hdd$counts to NA if they are zero
 newCounts<-hdd$counts
 newCounts[newCounts<(1e-05)]<-NA
 hdd$counts<-newCounts
 hdd
}


summary.plot.for.srac2lrac.single.ext1<-function(srac2lracTableList,lracDNAStringSet,lrac.id,srac.id,blastn6.aug.df,contigSummData,sracBinMembershipList,Xlim,Ylim)
{
 require(Biostrings)
 require(scales)
 
 srac2lracTableS<-srac2lracTableList[[lrac.id]][[srac.id]]
 if(nrow(srac2lracTableS)>0)
 {
  lrac.fasta<-lracDNAStringSet[lrac.id]
  #--turn this into "segment data"
  #--calculate the midpoints of each alignment
  xad<-rowMeans(cbind(srac2lracTableS[,"send"],srac2lracTableS[,"sstart"]))
 
  lo<-layout(matrix(c(1,1,2,2,3,4,5,6,7,8,9,10),4,3,byrow=F),height=c(0.8,rep(0.7,2),1),width=rep(1,4))
  
  binNum<-substring(srac.id,5)
  blastn6.aug.df.cut<-blastn6.aug.df[which(blastn6.aug.df$lrac==lrac.id),]
  binInd<-which(blastn6.aug.df.cut$sracBin==binNum)
  kappaVal<-blastn6.aug.df.cut[binInd,"kappa"]
  hd<-hist(blastn6.aug.df.cut$kappa,100,plot=F)
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count)*1.15),las=1,xlab=paste("Kappa statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(kappaVal,max(hd$count)*1.10,kappaVal,max(hd$count),length=0.1)
  text(kappaVal,max(hd$count)*1.15,srac.id)
  mtext(side=2,at=max(hd$count)*1.15,line=4,"A.",cex=1,las=1)

  par(mar=c(5.1,6.1,1.5,2.1))
  plot(contigSummData[,c("lgcov","gc")],pch=16,cex=0.75,col="lightgrey",las=1,xlab="Log coverage",ylab="GC",xlim=Xlim,ylim=Ylim)
  for(curbin in names(sracBinMembershipList))
  {
   Plot_ConvexHull(contigSummData[sracBinMembershipList[[curbin]],"lgcov"],contigSummData[sracBinMembershipList[[curbin]],"gc"],alpha("lightgrey",0.25),3)
  }
  points(contigSummData[sracBinMembershipList[[srac.id]],c("lgcov","gc")],pch=16,cex=0.75,col=1)
  legend(Xlim[1],Ylim[2],srac.id,ncol=1,pch=16,col=1,cex=1,bty="n")
  mtext(side=2,at=(Ylim[2]*1.01),line=4,"B.",cex=1,las=1)

  par(mar=c(0.1,6.1,1.5,2.1))
  #--reuse 'srac2lracTableS.segmentData' for each of the first 3 plots
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=srac2lracTableS[,"pident"],x1=srac2lracTableS[,"send"],y1=srac2lracTableS[,"pident"])
  plot(xad,rep(100,nrow(srac2lracTableS)),type="n",xaxt="n",las=1,ylab="PID",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=100,line=4,"C.",cex=1,las=1)
  
  par(mar=c(0.1,6.1,0.1,2.1))
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=srac2lracTableS[,"al2ql"],x1=srac2lracTableS[,"send"],y1=srac2lracTableS[,"al2ql"])
  plot(xad,rep(100,nrow(srac2lracTableS)),type="n",xaxt="n",las=1,ylab="al2ql",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1.12))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=1.1,line=4,"D.",cex=1,las=1)
  
  par(mar=c(0.1,6.1,0.1,2.1))
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=log10(srac2lracTableS[,"bitscore"]),x1=srac2lracTableS[,"send"],y1=log10(srac2lracTableS[,"bitscore"]))
  plot(xad,rep(6.5,nrow(srac2lracTableS)),type="n",xaxt="n",las=1,ylab="Bitscore (log10)",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(1.5,6.5))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=6.5,line=4,"E.",cex=1,las=1)

  #--next calculate GC-content of LRAC
  lrac.fasta.gc<-calc.gc.content.across.lrac(lrac.fasta,0.01)

  par(mar=c(5.1,6.1,0.1,2.1))
  plot(lrac.fasta.gc,las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="GC",pch=16,cex=0.90,col="lightgrey",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1))
  points(xad,srac2lracTableS[,"qgc"],col="black",pch=16)
  mtext(side=2,at=1,line=4,"F.",cex=1,las=1)
  legend(0,1,c(lrac.id,srac.id),ncol=2,bty="n",pch=16,col=c("lightgrey","black"))
  
  
  par(mar=c(0.1,6.1,1.5,2.1))
  numSracVal<-blastn6.aug.df.cut[binInd,"num.srac"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$num.srac,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab="",ylab="N",main="",xlim=c(0,1.05),xaxt="n")
  arrows(numSracVal,max(hd$count,na.rm=T)*1.10,numSracVal,max(hd$count,na.rm=T),length=0.05)
  text(numSracVal,max(hd$count,na.rm=T)*1.15,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"G.",cex=1,las=1)
  mtext(side=3,at=0.5,line=0,"num.srac",cex=0.75)

  par(mar=c(0.1,6.1,1.1,2.1))
  pidentVal<-blastn6.aug.df.cut[binInd,"pident"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$pident,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab="",ylab="N",main="",xlim=c(0,1.05),xaxt="n")
  arrows(pidentVal,max(hd$count,na.rm=T)*1.10,pidentVal,max(hd$count,na.rm=T),length=0.05)
  text(pidentVal,max(hd$count,na.rm=T)*1.15,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"H.",cex=1,las=1)
  mtext(side=3,at=0.5,line=0,"pident",cex=0.75)
  
  par(mar=c(0.1,6.1,1.1,2.1))
  al2qlVal<-blastn6.aug.df.cut[binInd,"al2ql"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$al2ql,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab="",ylab="N",main="",xlim=c(0,1.05),xaxt="n",cex.axis=0.9)
  arrows(al2qlVal,max(hd$count,na.rm=T)*1.10,al2qlVal,max(hd$count,na.rm=T),length=0.05)
  text(al2qlVal,max(hd$count,na.rm=T)*1.15,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"I.",cex=1,las=1)
  mtext(side=3,at=0.5,line=0,"al2ql",cex=0.75)
  
  par(mar=c(5.1,6.1,1.1,2.1))
  gappVal<-blastn6.aug.df.cut[binInd,"gapp"]
  hd<-hist(blastn6.aug.df.cut$gapp,100,plot=F)
  hd1<-set.zero.values.to.NA.in.hist.counts(hd)
  plot(hd1$mids,hd1$count,type="h",ylim=c(0,max(hd1$count,na.rm=T)*1.15),las=1,xlab=paste("Kappa statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1))
  arrows(gappVal,max(hd1$count,na.rm=T)*1.10,gappVal,max(hd1$count,na.rm=T),length=0.05)
  text(gappVal,max(hd1$count,na.rm=T)*1.15,srac.id,cex=0.85)
  mtext(side=2,at=max(hd1$count,na.rm=T)*1.15,line=4,"J.",cex=1,las=1)
  mtext(side=3,at=0.5,line=0,"gapp",cex=0.85)
 }
 else
 if(nrow(srac2lracTableS)<1)
 {
  print("No alignments are available for this pair!\n")
 }
 par(mfrow=c(1,1))
 #data.frame(srac2lracTableS,xad=xad,stringsAsFactors=F)
 return(hd)
}

#--updated version that displays coverage profiles
#--sr2lracCov,lr2lracCov are defined in /Users/rbhwilliams/Dropbox/nanopore-runs/run_2017_02_11/paper/g3.revision.190219/results/coverage

summary.plot.for.srac2lrac.single.ext2<-function(srac2lracTableList,lracDNAStringSet,lrac.id,srac.id,blastn6.aug.df,contigSummData,sracBinMembershipList,sr2lracCov,lr2lracCov,Xlim,Ylim)
{
 require(Biostrings)
 require(scales)
 
 srac2lracTableS<-srac2lracTableList[[lrac.id]][[srac.id]]
 if(nrow(srac2lracTableS)>0)
 {
  lrac.fasta<-lracDNAStringSet[lrac.id]
  #--turn this into "segment data"
  #--calculate the midpoints of each alignment
  xad<-rowMeans(cbind(srac2lracTableS[,"send"],srac2lracTableS[,"sstart"]))
 
  lo<-layout(matrix(c(1,1,2,2,3,3,4,5,5,6,6,7,8,9,9,10,10,11),6,3,byrow=F),height=c(0.25,0.06,0.19,0.19,0.06,0.29),width=rep(1,3))

  #--left column starts here

  par(mar=c(5.1,6.1,1.0,1.1))
  binNum<-substring(srac.id,5)
  blastn6.aug.df.cut<-blastn6.aug.df[which(blastn6.aug.df$lrac==lrac.id),]
  binInd<-which(blastn6.aug.df.cut$sracBin==binNum)
  kappaVal<-blastn6.aug.df.cut[binInd,"kappa"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$kappa,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("Kappa statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(kappaVal,max(hd$count,na.rm=T)*1.05,kappaVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(kappaVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"A.",cex=1,las=1)

  par(mar=c(5.1,6.1,1.0,1.1))
  plot(contigSummData[,c("lgcov","gc")],pch=16,cex=0.75,col="lightgrey",las=1,xlab="Log coverage",ylab="GC",xlim=Xlim,ylim=Ylim)
  for(curbin in names(sracBinMembershipList))
  {
   Plot_ConvexHull(contigSummData[sracBinMembershipList[[curbin]],"lgcov"],contigSummData[sracBinMembershipList[[curbin]],"gc"],alpha("lightgrey",0.25),3)
  }
  points(contigSummData[sracBinMembershipList[[srac.id]],c("lgcov","gc")],pch=16,cex=0.75,col=1)
  legend(Xlim[1],Ylim[2]*1.05,srac.id,ncol=1,pch=16,col=1,cex=1,bty="n")
  mtext(side=2,at=(Ylim[2]*1.01),line=4,"B.",cex=1,las=1)

  par(mar=c(5.1,6.1,1.0,1.1))
  plot(lr2lracCov,type="p",lty=1,ylim=c(0,(max(lr2lracCov[,2],sr2lracCov[,2])*1.1)),pch=3,cex=0.80,col="darkgrey",cex.lab=0.90,cex.axis=0.9,las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="Coverage")
  points(sr2lracCov,pch=3,cex=0.80,col="black",cex.lab=0.80,cex.axis=0.9,)
  abline(h=0,lty=3)
  legend(lr2lracCov[1,1],max(lr2lracCov[,2],sr2lracCov[,2])*1.15,c("SR","LR"),pch=3,cex=0.80,col=c("black","darkgrey"),ncol=2,bty="n")
  mtext(side=2,at=(max(lr2lracCov[,2],sr2lracCov[,2])*1.1),line=4,"C.",cex=1,las=1)
  
  #--middle column start here
  
  par(mar=c(4.1,6.1,1.0,1.1))
  #--reuse 'srac2lracTableS.segmentData' for each of the first 3 plots
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=srac2lracTableS[,"pident"],x1=srac2lracTableS[,"send"],y1=srac2lracTableS[,"pident"])
  plot(xad,rep(100,nrow(srac2lracTableS)),type="n",las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="PID",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=100,line=4,"D.",cex=1,las=1)
  
  par(mar=c(4.1,6.1,1.1,1.1))
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=srac2lracTableS[,"al2ql"],x1=srac2lracTableS[,"send"],y1=srac2lracTableS[,"al2ql"])
  plot(xad,rep(100,nrow(srac2lracTableS)),type="n",las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="al2ql",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1.12))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=1.1,line=4,"E.",cex=1,las=1)
  
  par(mar=c(4.1,6.1,1.1,1.1))
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=log10(srac2lracTableS[,"bitscore"]),x1=srac2lracTableS[,"send"],y1=log10(srac2lracTableS[,"bitscore"]))
  plot(xad,rep(6.5,nrow(srac2lracTableS)),type="n",las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="Bitscore (log10)",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(1.5,6.5))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=6.5,line=4,"F.",cex=1,las=1)

  #--next calculate GC-content of LRAC
  lrac.fasta.gc<-calc.gc.content.across.lrac(lrac.fasta,0.01)

  par(mar=c(5.1,6.1,1.1,1.1))
  plot(lrac.fasta.gc,las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="GC",pch=16,cex=0.90,col="lightgrey",cex.lab=0.90,cex.axis=0.9,xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1))
  points(xad,srac2lracTableS[,"qgc"],col="black",pch=16)
  mtext(side=2,at=1,line=4,"G.",cex=1,las=1)
  legend(0,1,c(lrac.id,srac.id),ncol=2,bty="n",pch=16,col=c("lightgrey","black"))
  
  #--right columns starts here
  
  par(mar=c(4.1,6.1,1.0,1.1))
  numSracVal<-blastn6.aug.df.cut[binInd,"num.srac"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$num.srac,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("num.srac statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(numSracVal,max(hd$count,na.rm=T)*1.05,numSracVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(numSracVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"H.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"num.srac",cex=0.75)

  par(mar=c(4.1,6.1,1.1,1.1))
  pidentVal<-blastn6.aug.df.cut[binInd,"pident"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$pident,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("pident statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(pidentVal,max(hd$count,na.rm=T)*1.05,pidentVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(pidentVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"I.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"pident",cex=0.75)
  
  par(mar=c(4.1,6.1,1.1,1.1))
  al2qlVal<-blastn6.aug.df.cut[binInd,"al2ql"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$al2ql,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("al2ql statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(al2qlVal,max(hd$count,na.rm=T)*1.05,al2qlVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(al2qlVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"J.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"al2ql",cex=0.75)
  
  par(mar=c(5.1,6.1,1.1,1.1))
  gappVal<-blastn6.aug.df.cut[binInd,"gapp"]
  hd<-hist(blastn6.aug.df.cut$gapp,100,plot=F)
  hd1<-set.zero.values.to.NA.in.hist.counts(hd)
  plot(hd1$mids,hd1$count,type="h",ylim=c(0,max(hd1$count,na.rm=T)*1.15),las=1,xlab=paste("gapp statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(gappVal,max(hd1$count,na.rm=T)*1.05,gappVal,max(hd1$count,na.rm=T)*0.90,length=0.05)
  text(gappVal,max(hd1$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd1$count,na.rm=T)*1.15,line=4,"K.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"gapp",cex=0.85)
 }
 else
 if(nrow(srac2lracTableS)<1)
 {
  print("No alignments are available for this pair!\n")
 }
 par(mfrow=c(1,1))
 #data.frame(srac2lracTableS,xad=xad,stringsAsFactors=F)
 return(hd)
}

summary.plot.for.srac2lrac.single.ext3<-function(srac2lracTableList,lracDNAStringSet,lrac.id,srac.id,blastn6.aug.df,contigSummData,sracBinMembershipList,sr2lracCov,lr2lracCov,Xlim,Ylim,covMax)
{
 require(Biostrings)
 require(scales)
 
 srac2lracTableS<-srac2lracTableList[[lrac.id]][[srac.id]]
 if(nrow(srac2lracTableS)>0)
 {
  lrac.fasta<-lracDNAStringSet[lrac.id]
  #--turn this into "segment data"
  #--calculate the midpoints of each alignment
  xad<-rowMeans(cbind(srac2lracTableS[,"send"],srac2lracTableS[,"sstart"]))
 
  lo<-layout(matrix(c(1,1,2,2,3,3,4,5,5,6,6,7,8,9,9,10,10,11),6,3,byrow=F),height=c(0.25,0.06,0.19,0.19,0.06,0.29),width=rep(1,3))

  #--left column starts here

  par(mar=c(5.1,6.1,1.0,1.1))
  binNum<-substring(srac.id,5)
  blastn6.aug.df.cut<-blastn6.aug.df[which(blastn6.aug.df$lrac==lrac.id),]
  binInd<-which(blastn6.aug.df.cut$sracBin==binNum)
  kappaVal<-blastn6.aug.df.cut[binInd,"kappa"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$kappa,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("Kappa statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(kappaVal,max(hd$count,na.rm=T)*1.05,kappaVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(kappaVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"A.",cex=1,las=1)

  par(mar=c(5.1,6.1,1.0,1.1))
  plot(contigSummData[,c("lgcov","gc")],pch=16,cex=0.75,col="lightgrey",las=1,xlab="Log coverage",ylab="GC",xlim=Xlim,ylim=Ylim*1.05)
  for(curbin in setdiff(names(sracBinMembershipList),srac.id))
  {
   Plot_ConvexHull(contigSummData[sracBinMembershipList[[curbin]],"lgcov"],contigSummData[sracBinMembershipList[[curbin]],"gc"],alpha("lightgrey",0.1),3)
  }
  Plot_ConvexHull(contigSummData[sracBinMembershipList[[srac.id]],"lgcov"],contigSummData[sracBinMembershipList[[srac.id]],"gc"],alpha("darkgrey",0.5),3)
  points(contigSummData[sracBinMembershipList[[srac.id]],c("lgcov","gc")],pch=16,cex=0.75,col=1)
  legend(Xlim[1],Ylim[2]*1.05,srac.id,ncol=1,pch=16,col=1,cex=0.85,bty="n")
  mtext(side=2,at=(Ylim[2]*1.01),line=4,"B.",cex=1,las=1)
 
  par(mar=c(5.1,6.1,1.0,1.1))
  covYlimMax<-pmin(covMax,max(lr2lracCov[,2],sr2lracCov[,2]))*1.1
  plot(lr2lracCov,type="p",lty=1,ylim=c(0,covYlimMax),pch=3,cex=0.80,col="darkgrey",cex.lab=0.90,las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="Coverage")
  points(sr2lracCov,pch=3,cex=0.80,col="black",cex.lab=0.80,cex.axis=0.9,)
  abline(h=0,lty=3)
  legend(lr2lracCov[1,1],covYlimMax,c("SR","LR"),pch=3,cex=0.85,col=c("black","darkgrey"),ncol=2,bty="n")
  mtext(side=2,at=(max(lr2lracCov[,2],sr2lracCov[,2])*1.1),line=4,"C.",cex=1,las=1)
  
  #--middle column start here
  
  par(mar=c(4.1,6.1,1.0,1.1))
  #--reuse 'srac2lracTableS.segmentData' for each of the first 3 plots
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=srac2lracTableS[,"pident"],x1=srac2lracTableS[,"send"],y1=srac2lracTableS[,"pident"])
  plot(xad,rep(100,nrow(srac2lracTableS)),type="n",las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="PID",xlim=c(0,nchar(lrac.fasta)),ylim=c(0,100))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=100,line=4,"D.",cex=1,las=1)
  
  par(mar=c(4.1,6.1,1.1,1.1))
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=srac2lracTableS[,"al2ql"],x1=srac2lracTableS[,"send"],y1=srac2lracTableS[,"al2ql"])
  plot(xad,rep(100,nrow(srac2lracTableS)),type="n",las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="al2ql",xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1.12))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=1.1,line=4,"E.",cex=1,las=1)
  
  par(mar=c(4.1,6.1,1.1,1.1))
  srac2lracTableS.segmentData<-data.frame(x0=srac2lracTableS[,"sstart"],y0=log10(srac2lracTableS[,"bitscore"]),x1=srac2lracTableS[,"send"],y1=log10(srac2lracTableS[,"bitscore"]))
  plot(xad,rep(6.5,nrow(srac2lracTableS)),type="n",las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="Bitscore (log10)",xlim=c(0,nchar(lrac.fasta)),ylim=c(1.5,6.5))
  plot.segments(srac2lracTableS.segmentData,col=1,lwd=2)
  mtext(side=2,at=6.5,line=4,"F.",cex=1,las=1)

  #--next calculate GC-content of LRAC
  lrac.fasta.gc<-calc.gc.content.across.lrac(lrac.fasta,0.01)

  par(mar=c(5.1,6.1,1.1,1.1))
  plot(lrac.fasta.gc,las=1,xlab=paste("Physical location on LR-chr",lrac.id,"(bp)",sep=" "),ylab="GC",pch=16,cex=0.90,col="lightgrey",xlim=c(0,nchar(lrac.fasta)),ylim=c(0,1))
  points(xad,srac2lracTableS[,"qgc"],col="black",pch=16)
  mtext(side=2,at=1,line=4,"G.",cex=1,las=1)
  legend(0,1,c(lrac.id,srac.id),ncol=2,bty="n",pch=16,col=c("lightgrey","black"),cex=0.85)
  
  #--right columns starts here
  
  par(mar=c(4.1,6.1,1.0,1.1))
  numSracVal<-blastn6.aug.df.cut[binInd,"num.srac"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$num.srac,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("num.srac statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(numSracVal,max(hd$count,na.rm=T)*1.05,numSracVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(numSracVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"H.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"num.srac",cex=0.75)

  par(mar=c(4.1,6.1,1.1,1.1))
  pidentVal<-blastn6.aug.df.cut[binInd,"pident"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$pident,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("pident statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(pidentVal,max(hd$count,na.rm=T)*1.05,pidentVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(pidentVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"I.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"pident",cex=0.75)
  
  par(mar=c(4.1,6.1,1.1,1.1))
  al2qlVal<-blastn6.aug.df.cut[binInd,"al2ql"]
  hd<-set.zero.values.to.NA.in.hist.counts(hist(blastn6.aug.df.cut$al2ql,100,plot=F))
  plot(hd$mids,hd$count,type="h",ylim=c(0,max(hd$count,na.rm=T)*1.15),las=1,xlab=paste("al2ql statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(al2qlVal,max(hd$count,na.rm=T)*1.05,al2qlVal,max(hd$count,na.rm=T)*0.90,length=0.05)
  text(al2qlVal,max(hd$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd$count,na.rm=T)*1.15,line=4,"J.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"al2ql",cex=0.75)
  
  par(mar=c(5.1,6.1,1.1,1.1))
  gappVal<-blastn6.aug.df.cut[binInd,"gapp"]
  hd<-hist(blastn6.aug.df.cut$gapp,100,plot=F)
  hd1<-set.zero.values.to.NA.in.hist.counts(hd)
  plot(hd1$mids,hd1$count,type="h",ylim=c(0,max(hd1$count,na.rm=T)*1.15),las=1,xlab=paste("gapp statistics for",lrac.id,sep=" "),ylab="N",main="",xlim=c(0,1.05))
  arrows(gappVal,max(hd1$count,na.rm=T)*1.05,gappVal,max(hd1$count,na.rm=T)*0.90,length=0.05)
  text(gappVal,max(hd1$count,na.rm=T)*1.10,srac.id,cex=0.85)
  mtext(side=2,at=max(hd1$count,na.rm=T)*1.15,line=4,"K.",cex=1,las=1)
  #mtext(side=3,at=0.5,line=0,"gapp",cex=0.85)
 }
 else
 if(nrow(srac2lracTableS)<1)
 {
  print("No alignments are available for this pair!\n")
 }
 par(mfrow=c(1,1))
 #data.frame(srac2lracTableS,xad=xad,stringsAsFactors=F)
 return(hd)
}

