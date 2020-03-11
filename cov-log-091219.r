#--plotting up some the coverage data

#--take a look at PAO1-tig*1
#lr.pao1.tig00000001<-read.table(file="minimap2_LR_LRAC.corr/PAO1-tig00000001_perbasecov",sep="\t",header=F)
#sr.pao1.tig00000001<-read.table(file="minimap2_SR_LRAC.corr/PAO1-tig00000001_perbasecov",sep="\t",header=F)

#--median LR coverage
#> median(lr.pao1.tig00000001[,3])
#[1] 400

#--number of zeros?
#> length(which(lr.pao1.tig00000001[,3]<1))
#[1] 0

#--median SR coverage
#> median(sr.pao1.tig00000001[,3])
#[1] 757

#--and number of zeros
#> length(which(sr.pao1.tig00000001[,3]<1))
#[1] 858

#--modified on 18 Dec 2019; copy this to a separete file to send to Mindia
count.zero.segments<-function(x,thrsh=1,rank=T)
{
 #--'x' is coverage vector,
 #--'thrsh' is the threshold around which the signal is hard threshold to 0 or 1
 #--if 'rank' is TRUE, the results will be sorted by the inteval in decreasing order
 y<-rep(0,length(x))
 y[x>=thrsh]<-1
 #--extend to avoid end-zeros
 y<-c(1,y,1)
 #--and take differences
 yd<-diff(y)
 #--define intervals and duration of runs of zeros
 startInd<-which(yd<0)
 endInd<-which(yd>0)-1
 duration<-endInd-startInd+1
 res<-data.frame(start=startInd,end=endInd,interval=duration,stringsAsFactors=F)
 #--if rank==T, sort rows by interval in decreasing order
 if(rank==T)
 {
  res<-res[order(res$interval,decreasing=T),]
 }
 return(res)
}

#--use this as a test case
#testseq<-c(0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0,1)

#> testseq
#[1] 0 1 1 1 0 0 0 0 0 0 1 1 1 0 0 1 1 1 0 0 1

#> count.zero.segments(testseq,1)
#[1] 1 6 2 2

#--and then on the short read data
#> count.zero.segments(sr.pao1.tig00000001[,3],1)
#[1]   1   1  60   2   5 116   1   4   2 186  17  73  29   1   1   1  42  38   1
#[20]  42   3 155  77

#--so we can extract statisics such as the maximum run of zero-coverage etc.

#--Rohan finished editing here on 10 December 2019--

#--check the place that Mindi mentioned is possible mis-assembled?
#plot(lr.pao1.tig00000001[1415000:1420000,2:3],ylim=c(0,1000),pch=16,cex=0.25,las=1,col="darkgrey",xlab="Position on PAO1-tig00000001 (bp)",ylab="Coverage")
#points(sr.pao1.tig00000001[1415000:1420000,2:3],pch=16,cex=0.25,col="black")
#abline(h=0,lty=2)
#legend(1415000,1000,c("LR","SR"),col=c("black","darkgrey"),pch=16,cex=0.80)

#--added on 12 December 2019--
#--process coverage profiles so we can use them in srac2lrac panel plots

reduce.coverage.data<-function(X,d=0.01)
{
 #--etimate window size k
 w<-min(diff(round(quantile(1:length(X),probs=seq(0,1,d)),0)))
 Xrmed<-runmed(X,w)
 
 #--define summary points
 p<-seq(w/2,length(X),by=w)
 
 #--and output
 res<-data.frame(pos=p,cov=Xrmed[p])
 return(res)
}

#--repeat for all 21 genomes
#lr2lracFiles<-dir("minimap2_LR_LRAC.corr/",patt="perbasecov")
#lr2lracCovList<-list()
#for(curfile in lr2lracFiles)
#{
# curpath<-paste("minimap2_LR_LRAC.corr/",curfile,sep="")
# curdata<-read.table(curpath,sep="\t",header=F)
# curtag<-strsplit(curfile,"_")[[1]][1]
# print(curtag)
# lr2lracCovList[[curtag]]<-reduce.coverage.data(curdata[,3],0.0025)
#}

#--repeat for SR data
#sr2lracFiles<-dir("minimap2_SR_LRAC.corr/",patt="perbasecov")
#sr2lracCovList<-list()
#for(curfile in sr2lracFiles)
#{
# curpath<-paste("minimap2_SR_LRAC.corr/",curfile,sep="")
# curdata<-read.table(curpath,sep="\t",header=F)
# curtag<-strsplit(curfile,"_")[[1]][1]
# print(curtag)
# sr2lracCovList[[curtag]]<-reduce.coverage.data(curdata[,3],0.0025)
#}

#save(lr2lracCovList,sr2lracCovList,file="reducedCoverageProfiles.RData")

#--added on 19 Dec 2019
#--update this function to detect both runs of both low and high values

source(file="find.outlier.segments.r")

#testseq1<-c(0,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,0,0,1,1,1,0,0,1)

#> get.outlier.segments(testseq1,1,"below",F)
#start end interval
#1     1   1        1
#2     5  10        6
#3    19  20        2
#4    24  25        2

#--added on February 26 2020
#--Krithika has computed short and long read mapping for "artefact" case
#./minimap2_LR_LRAC.corr/artefact/PAO3A-tig00000001_perbasecov
#./minimap2_SR_LRAC.corr/artefact/PAO3A-tig00000001_perbasecov
#--so read in and process as above

#lr2lracArtefactList<-list()
#artefactExample.lr2lracCov<-read.table("./minimap2_LR_LRAC.corr/artefact/PAO3A-tig00000001_perbasecov",sep="\t",header=F)
#lr2lracArtefactList[["PAO3A-tig00000001"]]<-reduce.coverage.data(artefactExample.lr2lracCov[,3],0.0025)

#sr2lracArtefactList<-list()
#artefactExample.sr2lracCov<-read.table("./minimap2_SR_LRAC.corr/artefact/PAO3A-tig00000001_perbasecov",sep="\t",header=F)
#sr2lracArtefactList[["PAO3A-tig00000001"]]<-reduce.coverage.data(artefactExample.sr2lracCov[,3],0.0025)

#save(lr2lracArtefactList,sr2lracArtefactList,file="reducedCoverageProfilesArtefacts.RData")

#--added on 27 February 2020
#--Krithika has forwarded the per base coverage from new SR data from PAO4, so update 'sr2lracCovList'
pao4.perbasecov.files<-list.files(path="minimap2_SR_LRAC.corr/",patt="PAO4")
for(curfile in pao4.perbasecov.files)
{
 curpath<-paste("minimap2_SR_LRAC.corr/",curfile,sep="")
 curdata<-read.table(curpath,sep="\t",header=F)
 curtag<-strsplit(curfile,"_")[[1]][1]
 print(curtag)
 sr2lracCovList[[curtag]]<-reduce.coverage.data(curdata[,3],0.0025)
}

#--and resave
save(lr2lracCovList,sr2lracCovList,file="reducedCoverageProfiles.RData")
