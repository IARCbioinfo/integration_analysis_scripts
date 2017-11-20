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
require(gplots)
require(lattice)
# read data
Ddata = NULL
Cdata = NULL
Rdata = NULL
Mdata = NULL


dtl = vector("list",4)
ntypes = 0
type = c()
row.order = c()
plot.chr  = c()
sparse    = c()
cap       = c()
if(!is.null(opt$DNA)){
  ntypes = ntypes +1
  Ddata = read.table(opt$DNA,h=T,row.names = 1)
  dtl[[ntypes]] = t(Ddata)
  type = c(type,"binomial")
  row.order = c(row.order,T)
  plot.chr  = c(plot.chr,F)
  sparse    = c(sparse,T)
  cap       = c(cap,T)
}
if(!is.null(opt$CNV)){
  ntypes = ntypes +1
  dtl[[ntypes]] = t(Cdata)
  type = c(type,"gaussian")
  row.order = c(row.order,T)
  plot.chr  = c(plot.chr,T)
  sparse    = c(sparse,T)
  cap       = c(cap,T)
}
if(!is.null(opt$RNA)){
  ntypes = ntypes +1
  Rdata = read.table(opt$RNA,h=T,row.names = 1,check.names = T)
  dtl[[ntypes]] = t(Rdata)
  type = c(type,"gaussian")
  row.order = c(row.order,T)
  plot.chr  = c(plot.chr,F)
  sparse    = c(sparse,T)
  cap       = c(cap,T)
}
if(!is.null(opt$methylation)){
  ntypes = ntypes +1
  Mdata = read.table(opt$methylation,h=T,row.names = 1)
  dtl[[ntypes]] = t(Mdata)
  type = c(type,"gaussian")
  row.order = c(row.order,T)
  plot.chr  = c(plot.chr,F)
  sparse    = c(sparse,T)
  cap       = c(cap,T)
}

print(dtl)
print(type)

# tune lambda
cv.fitl = vector("list",Kmax-Kmin+1)
for(k in (Kmin:Kmax-1) ){
  cv.fitl[[k]] = tune.iClusterPlus(cpus=cores,dt1=dtl[[1]],dt2=dtl[[2]],dt3=dtl[[3]],dt4=dtl[[4]], type=type, n.lambda=NULL,K=k,maxiter=20)
  cv.fittmp = cv.fitl[[k]]
  save(cv.fittmp, file=paste(opt$out,"cv.fit.K",k+1,".Rdata",sep=""))
}
warnings()
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
clusters2 = clusters
k = Kmin
while(k<Kmax){
  print(k)
  clusters2 = clusters
  for(i in 1:max(clusters[,k-1]) ){
    print(i)
    tatmp = table(clusters[clusters[,k-1]==i,k])
    newid = as.numeric(names(tatmp)[which.max(tatmp)])
    clusters2[clusters[,k]==newid,k] = i
    clusters2[clusters[,k]==i,k] = newid
    print(clusters2)
    clusters = clusters2
  }
  k = k+1
}

svg(paste(opt$out,"iClusters.svg",sep=""),h=6,w=8*2)
par(family="Times",mfrow=c(1,2))
plot(Kmin:Kmax,devRatMinBIC,type="b",xlab=expression(italic(K)),ylab=expression("Pseudo "~italic(R^2)),las=1)
image(1:nrow(clusters),Kmin:Kmax,clusters,col=rainbow(Kmax),axes=F,xlab="",ylab=expression(italic(K)),las=2 )
axis(1,at = 1:nrow(clusters),labels = rownames(dtl[[1]]) ,las=2)
axis(2,at = Kmin:Kmax,las=1)
dev.off()

# plot heatmap and write features
features = lapply(dtl,colnames)
for(k in (Kmin:Kmax-1) ){
  print(k)
  best.cluster=clusters[,k]
  print(best.cluster)
  best.fit=cv.fitl[[k]]$fit[[which.min(BIC[,k])]]
  print(best.fit)
  col.scheme = alist()
  col.scheme[[1]] = bluered(256)
  col.scheme[[2]] = bluered(256)
  col.scheme[[3]] = bluered(256)
  col.scheme[[4]] = bluered(256)

  png(paste(opt$out,"Heatmap_K",k+1,".png",sep=""),h=4*200,w=4*300,bg=0)
  plotHeatmap(fit=best.fit,datasets=dtl[1:ntypes], type=type, col.scheme = col.scheme[1:ntypes], row.order=row.order,plot.chr=plot.chr,sparse=sparse,cap=cap)
  dev.off()
  
  # get best features
  sigfeatures=alist()
  idsigfeatures=alist()
  print("yop")
  print(sum(!sapply(features,is.null)))
  for(i in 1:sum(!sapply(features,is.null)) ){
    print(i)
    rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
    upper=quantile(rowsum,prob=0.75)
    sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
    idsigfeatures[[i]]=which(rowsum>upper)
    write.table(sigfeatures[[i]],file = paste(opt$out,"_K",k+1,"_features",i,".txt",sep=""),col.names = F,quote = F)
  }
  #names(sigfeatures)=c("expression","methylation")
}



