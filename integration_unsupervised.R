##########################################################################
### Script to perform integrative Clustering of several omics datasets ###
##########################################################################

# get options
library("optparse")

option_list = list(
  make_option(c("-d", "--DNA"), type="character", default=NULL, help="file with mutation data [default= %default]", metavar="character"),
  make_option(c("-r", "--RNA"), type="character", default=NULL, help="file with expression data [default= %default]", metavar="character"),
  make_option(c("-k", "--Kmin"), type="numeric", default=2, help="minimum number of clusters [default= %default]", metavar="numeric"),
  make_option(c("-K", "--Kmax"), type="numeric", default=6, help="maximum number of clusters [default= %default]", metavar="numeric"),
  make_option(c("-m", "--methylation"), type="character", default=NULL, help="file with methylation data [default= %default]", metavar="character"),
  make_option(c("-C", "--CNV"), type="character", default=NULL, help="file with copy number data [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, help="number of cores for statistical computation [default= %default]", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Kmin  = as.numeric(opt$Kmin)
Kmax  = as.numeric(opt$Kmax)
cores = as.numeric(opt$cores)

require(iClusterPlus)

# read data
if(!is.null(opt$DNA)) Ddata = read.table(opt$DNA,h=T)
if(!is.null(opt$RNA)) Rdata = read.table(opt$RNA,h=T)
if(!is.null(opt$methylation)) Mdata = read.table(opt$Methylation,h=T)
if(!is.null(opt$CNV)) Cdata = read.table(opt$CNV,h=T)

# tune lambda
cv.fitl = vector("list",5)
for(k in (Kmin:Kmax-1) ){
  cv.fitl[[k]] = tune.iClusterPlus(cpus=cores,dt1=t(Ddata),dt2=t(Cdata),dt3=t(Rdata),dt4=t(Mdata), type=c("binomial","gaussian","gaussian","gaussian"), n.lambda=NULL,K=k,maxiter=20)
  cv.fittmp = cv.fitl[[k]]
  save(cv.fittmp, file=paste(opt$out,"cv.fit.k",k,".Rdata",sep=""))
}

nLambda = nrow(cv.fitl[[1]]$lambda)
nK = length(cv.fitl)
BIC = getBIC(cv.fitl)
devR = getDevR(cv.fitl)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

#best.fitl=sapply( Kmin:Kmax-1 , function(k) cv.fitl[[k]]$fit[[which.min(BIC[,k])]] )
clusters=getClusters(cv.fitl)
colnames(clusters)=paste("K=",Kmin:Kmax,sep="")


svg(paste(out,"iClusters.svg",sep=""),h=6,w=8*2)
par(family="Times",las=1,mfrow=c(1,2))
plot(Kmin:Kmax,devRatMinBIC,type="b",xlab="K",ylab="Pseudo R2")
image(cbind(clusters),col=prettycolors,axes=T,xlab="",ylab="K")
#axis(1,at = 1:50,labels=rep("",50) )
#axis(2,at = 2:6,las=1)
dev.off()

# plot heatmap
for(k in (Kmin:Kmax-1) ){
  best.cluster=clusters[,k]
  best.fit=cv.fitl[[k]]$fit[[which.min(BIC[,k])]]
  col.scheme = alist()
  col.scheme[[1]] = bluered(256)
  col.scheme[[2]] = bluered(256)
  col.scheme[[3]] = bluered(256)
  col.scheme[[4]] = bluered(256)

  png(paste("Heatmap_K",k+1,".png",sep=""),h=4*200,w=4*300,bg=0)
  plotHeatmap(fit=best.fit,datasets=list(mutations=t(Ddata),CNV=t(Cdata),RNA=t(Rdata),Methylation=t(Mdata)), type=c("binomial","gaussian","gaussian","gaussian"), col.scheme = col.scheme, row.order=c(T,T,T,T),plot.chr=c(F,T,F,F),sparse=c(T,T,T,T),cap=c(T,T,T,T))
  dev.off()
}

# get best features
features = alist()
features[[1]] = rownames(Ddata)
features[[2]] = rownames(Cdata)
features[[3]] = rownames(Rdata)
features[[4]] = rownames(Mdata)
sigfeatures=alist()
idsigfeatures=alist()
for(i in 1:4){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
  idsigfeatures[[i]]=which(rowsum>upper)
}
names(sigfeatures)=c("expression","methylation")

write.table(sigfeatures,file = paste(out,"_features",sep=""))

