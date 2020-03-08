#--srac2lrac: functions for computing concordnace statistic between short read MAG and LR-chr sequences
#--this R code by Rohan Williams (rbhwilliams@gmail.com)
#--original developed as 'augment.srac2lrac.tables.v3' in ~/Documents/code/lrac2sracR/
#--renamed to 'augment.srac2lrac.tables' here
#--the fuction 'calculate.gaps.simple' was corrected by Rohan on 15 August 2019 to avoid the bug when using unicycler LRAC named as integers.
#--this file updated 16 August 2019; all function names remain unchanged.
#--updated to include 'augment.kappa.table' on 24 October 2019

augment.srac2lrac.tables<-function(srac2lracTable,lracSeq,sracSummTable)
{
 #--add some columns containing SRAC and LRAC inforamtion to a BLASTN table
 #--'srac2lracTable' is a data.frame containing the output of a BLASTN analysis (column should be named)
 #--'lracSeq' is a Biostrings DNAStringSet containing the LRAC sequences (subjects)
 #--'sracSummTable' is an output from RKXM::process.spades.contigs
 res<-srac2lracTable

 #--obtain the query (srac) lengtha
 qlen<-sracSummTable[srac2lracTable$qseqid,"len"]

 #--calculate quotient of alignment to query (srac) length
 al2ql<-srac2lracTable$length/qlen

 #--and add summary data from binning variables
 qlgcov<-log10(sracSummTable[srac2lracTable$qseqid,"cov"])
 qgc<-sracSummTable[srac2lracTable$qseqid,"gc"]

 res<-data.frame(res,qlen=qlen,al2ql=al2ql,qlgcov=qlgcov,qgc=qgc,stringsAsFactors=F)
 res
}

#--filter augmented BLASTN results for best hit
#--original developed as 'filter.hits.2' in ~/Documents/code/lrac2sracR/
filter.hits<-function(srac2lracTable)
{
 #--this function filters for best hit within each query-subject combination
 combTags<-apply(cbind(srac2lracTable[,1:2]),1,paste,collapse="-")
 ucombTags<-unique(combTags)
 res<-srac2lracTable[match(ucombTags,combTags),]
 res
}

#--this function generates the aggregated results that are used to explore associations between LRAC and SRAC (as grouped by bin/contig set)
build.results.cube.extended<-function(srac2lracTable,lrac.fasta,binMembershipList)
{
 #--make the output
 res<-list(sracSets=NULL,sracCaptureStats=NULL,cube=NULL,data=NULL)

 #--generate contig sets
 sracWithHits<-unique(srac2lracTable$qseqid)
 lracWithHits<-unique(srac2lracTable$sseqid)
 blastnVariables<-c("num.hits","num.srac","pident","length","mismatch","gapopen","evalue","bitscore","al2ql","gapp")
 
 res$sracSets<-lapply(binMembershipList,FUN=intersect,y=sracWithHits)
 res$data<-list()
 
 #--and contig set statistics relative to totals
 res$sracCaptureStats<-cbind(sapply(res$sracSets,length),sapply(binMembershipList,length))
 colnames(res$sracCaptureStats)<-c("with.hits","total")
 rownames(res$sracCaptureStats)<-names(res$sracSets)

 #--generate main set of results
 res$cube<-array(NA,dim=c(length(lracWithHits),length(res$sracSets),length(blastnVariables)),dimnames=list(lracWithHits,names(res$sracSets),blastnVariables))

 for(cur.lrac in dimnames(res$cube)[[1]])
 {
  print(cur.lrac)
  cur.srac2lracTable<-subset(srac2lracTable,sseqid==cur.lrac)
  if(nrow(cur.srac2lracTable)>1)
  {
   for(cur.sracSet in dimnames(res$cube)[[2]])
   {
    curind<-which(cur.srac2lracTable$qseqid%in%res$sracSets[[cur.sracSet]])
    if(length(curind)>1)
    {
     #--calculate blastn summary statistics
     cur.blastnResults<-cur.srac2lracTable[curind,]
     cur.blastnStats<-apply(cur.blastnResults[,c("pident","length","mismatch","gapopen","evalue","bitscore","al2ql")],2,median)
     res$data[[cur.lrac]][[cur.sracSet]]<-cur.blastnResults
     #--and basic srac statistics
     cur.num.hits<-length(curind)
     cur.num.srac<-length(unique(cur.srac2lracTable[curind,"qseqid"]))
     #--and gap statistics
     cur.gapp<-calculate.gaps.simple(cur.blastnResults,lrac.fasta)
     #--construct output
     cur.res<-c(cur.num.hits,cur.num.srac,cur.blastnStats,cur.gapp)
     res$cube[cur.lrac,cur.sracSet,]<-cur.res
    }
   }
  }
 }
 res
}

#--this function is used by 'build.results.cube.extended'
calculate.gaps.simple<-function(srac2lracTable,lrac.fasta)
{
 #--uses a simple direct calculation of proportion of an LRAC sequence that has SRAC alignments
 #--the alignment data is in 'srac2lracTable' (and MUST contain only results from a single LRAC sequence) and 'lrac.fasta' are the LRAC sequences captured as an DNAStringSet
 lracid<-as.character(unique(srac2lracTable[,"sseqid"]))#<--error was fixed by forcing lracid toa character not an integer
 stopifnot(length(lracid)==1)
 {
  posvec<-rep_len(x=0,length.out=nchar(lrac.fasta[lracid])) #<--error is coming in here, why? lracid is obtained from srac2lracTable, so it must be a problem in that input argument!
  #--interpolate alignments
  for(currow in (1:nrow(srac2lracTable)))
  {
   posvec[srac2lracTable[currow,"sstart"]:srac2lracTable[currow,"send"]]<-(posvec[srac2lracTable[currow,"sstart"]:srac2lracTable[currow,"send"]]+1)
  }
  res<-length(which(posvec>1e-15))/length(posvec)
 }
 return(res)
}

#--write a function to extract a simple summary statistic that seeks to identify {srac2lrac}|bin,lrac} that are full length, srac-bin complete, high pident and high alignment to query length ratio
cut.cube.for.kappa<-function(cube)
{
 require(abind)
 numSracNorm<-sweep(cube$cube[,rownames(cube$sracCaptureStats),"num.srac"],2,cube$sracCaptureStats[,"total"],"/")
 pidentNorm<-cube$cube[,rownames(cube$sracCaptureStats),"pident"]/100
 res<-abind(numSracNorm,pidentNorm,cube$cube[,rownames(cube$sracCaptureStats),"al2ql"],cube$cube[,rownames(cube$sracCaptureStats),"gapp"],along=3,new.names=c("num.srac","pident","al2ql","gapp"))
 res
}

calculate.kappa<-function(cube)
{
 res<-matrix(NA,nrow=dim(cube)[1],ncol=dim(cube)[2])
 rownames(res)<-dimnames(cube)[[1]]
 colnames(res)<-dimnames(cube)[[2]]
 for(curr in rownames(res))
 {
  print(curr)
  for(curc in colnames(res))
  {
   res[curr,curc]<-mean(cube[curr,curc,],na.rm=T)
  }
 }
 res
}

process.kappa.results<-function(mat)
{
 res=NULL
 curColNewName=NULL
 resName=NULL
 lracName=NULL
 sracBinName=NULL
 for(curcol in (1:ncol(mat)))
 {
  res<-c(res,as.vector(mat[,curcol]))
  curColNewName<-paste(rownames(mat),strsplit(colnames(mat)[curcol],".",fixed=T)[[1]][2],sep="-")
  resName<-c(resName,curColNewName)
  lracName<-c(lracName,rownames(mat))
  sracBinName<-c(sracBinName,rep(strsplit(colnames(mat)[curcol],".",fixed=T)[[1]][2],nrow(mat)))
 }
 res<-data.frame(tag=resName,lrac=lracName,sracBin=sracBinName,kappa=res,stringsAsFactors=F)
 res<-res[!is.na(res$kappa),]
 res<-res[order(res$kappa,decreasing=T),]
 rownames(res)<-res$tag
 return(res)
}

#--add a new function to that contructs a data cube indexed as lrac-by-srac(contigs)-by-blastn.statistics
build.contig.level.cube<-function(srac2lracTable)
{
 qidu<-unique(srac2lracTable$qseqid)
 sidu<-unique(srac2lracTable$sseqid)
 bstats.colInd<-c(3:6,11:14)
 bstats<-colnames(blastn2)[bstats.colInd]
 res<-array(NA,dim=c(length(sidu),length(qidu),length(bstats)),dimnames=list(lrac=sidu,srac=qidu,bnstat=bstats))
 
 #--loop and fill
 for(currow in (1:nrow(srac2lracTable)))
 {
  cat(currow)
  curq<-srac2lracTable$qseqid[currow]
  curs<-srac2lracTable$sseqid[currow]
  res[curs,curq,]<-as.numeric(srac2lracTable[currow,bstats.colInd])
 }
 res
}

#--added on 24 October 2019
#--add the '4' output to '6'
augment.kappa.table<-function(kappa.df,kappaStatsCube)
{
 #--augment the output of 'process.kappa.results' (kappa.df) to include the component-statistics of kappa generated by 'cut.cube.for.kappa' (kappaStatsCube)
 res=NULL
 augout=NULL
 for(currow in (1:nrow(kappa.df)))
 {
  curlrac<-kappa.df$lrac[currow]
  cursracBin<-paste("bin",kappa.df$sracBin[currow],sep=".")
  curres<-kappaStatsCube[curlrac,cursracBin,]
  augout<-rbind(augout,curres)
  print(currow)
 }
 res<-data.frame(kappa.df,augout)
 res
}
