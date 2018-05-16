#!/usr/bin/env Rscript

####################################################
### Script to perform Clustering of RNA seq data ###
####################################################

# get options
library("optparse")

option_list = list(
  make_option(c("-g", "--group_file"), type="character", default=NULL, help="file with sample groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="MOFA_out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-s", "--suffix"), type="character", default="count.txt", help="Suffix for output names [default= %default]", metavar="character"),
  make_option(c("-R", "--expression"), type="character", default=NULL, help="Expression dataset [default= %default]", metavar="character"),
  make_option(c("-M", "--methylation"), type="character", default=NULL, help="Methylation dataset [default= %default]", metavar="character"),
  make_option(c("-S", "--mutation"), type="character", default=NULL, help="Mutations dataset [default= %default]", metavar="character"),
  make_option(c("-r", "--robustness"), action="store_true", default=TRUE, help="Perform robustness analysis [default= %default]"),
  make_option(c("-i", "--maxiter"), type="numeric", default=10000, help="Maximum number of iterations [default= %default]"),
  make_option(c("-n", "--nconsreps"), type="numeric", default=100, help="Number of iterations for consensus clustering [default= %default]"),
  make_option(c("-m", "--nreps"), type="numeric", default=10, help="Maximum number of iterations for Latemt Factor robustness analysis [default= %default]")
); 

require(MOFAtools)
## add something to specify class of group file columns

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read data
data = list()
if(!(is.null(opt$mutation)))    data$Mutation = read.table(opt$mutation,h=T,row.names = 1)
if(!(is.null(opt$expression)))  data$RNA      = read.table(opt$expression,h=T,row.names = 1)
if(!(is.null(opt$methylation))) data$Methyl   = read.table(opt$methylation,h=T,row.names = 1)

if(!is.null(opt$robustness) ){
  print("Perform robustness analysis")
  nrep = as.numeric(opt$nconsreps)
  pItem = 0.8
  for(i in 1:nrep){
    print(i)
    print(c("outfile",paste(opt$out,"/MOFA",opt$suffix,"_sub",i,".hdf5",sep="")))
    rm( MOFAobjecttmp )
    MOFAobjecttmp <- createMOFAobject(data)
    nsamp = MOFAobjecttmp@Dimensions$N
    samptmp = sort(sample(1:nsamp,size = round(pItem*nsamp),replace=F) )
    MOFAobjecttmp@Dimensions$N = round(pItem*nsamp)
    for(j in 1:length(MOFAobjecttmp@TrainData) ) MOFAobjecttmp@TrainData[[j]] = MOFAobjecttmp@TrainData[[j]][,samptmp]
    
    DirOptions <- list( "dataDir" = tempdir(), "outFile" = paste(opt$out,"/MOFA",opt$suffix,"_sub",i,".hdf5",sep="") )
    ModelOptions <- getDefaultModelOptions(MOFAobjecttmp)
    ModelOptions$likelihood[names(ModelOptions$likelihood)=="Mutation"] = "bernoulli"
    ModelOptions$learnIntercept = F
    TrainOptions <- getDefaultTrainOptions()
    TrainOptions$maxiter = as.numeric(opt$maxiter)
    DataOptions <- getDefaultDataOptions()
    DataOptions$centerFeatures = T
    print(DirOptions)
    print(ModelOptions)
    print(TrainOptions)
    print(DataOptions)
    MOFAobjecttmp <- prepareMOFA(MOFAobjecttmp, DirOptions = DirOptions,ModelOptions = ModelOptions,TrainOptions = TrainOptions,DataOptions = DataOptions)
  
    # run
    MOFAobjecttmp <- runMOFA(MOFAobjecttmp, DirOptions)
  }
}

#perform MOFA run
print("Perform multiple runs")
MOFAfinal=c()
bestELBO = -Inf
bestmodel = 1
nrep = as.numeric(opt$nreps)
for(i in 1:nrep){
    print(i)
    print(c("outfile",paste(opt$out,"/MOFA",opt$suffix,"_run",i,".hdf5",sep="")) )
    MOFAobjecttmp <- createMOFAobject(data)
    DirOptions <- list( "dataDir" = tempdir(), "outFile" = paste(opt$out,"/MOFA",opt$suffix,"_run",i,".hdf5",sep="") )
    ModelOptions <- getDefaultModelOptions(MOFAobjecttmp)
    ModelOptions$likelihood[names(ModelOptions$likelihood)=="Mutation"] = "bernoulli"
    ModelOptions$learnIntercept = F
    TrainOptions <- getDefaultTrainOptions()
    TrainOptions$maxiter = as.numeric(opt$maxiter)
    DataOptions <- getDefaultDataOptions()
    DataOptions$centerFeatures = T
    print(DirOptions)
    print(ModelOptions)
    print(TrainOptions)
    print(DataOptions)
    MOFAobjecttmp <- prepareMOFA(MOFAobjecttmp, DirOptions = DirOptions,ModelOptions = ModelOptions,TrainOptions = TrainOptions,DataOptions = DataOptions)
    
    # run
    MOFAobjecttmp <- runMOFA(MOFAobjecttmp, DirOptions)
    
    if( rev(MOFAobjecttmp@TrainStats$elbo)[1]> bestELBO){
      bestELBO = rev(MOFAobjecttmp@TrainStats$elbo)[1]
      MOFAfinal = MOFAobjecttmp
      bestmodel = i
    }
}
MOFAobject = MOFAfinal
print(paste("Best model:",bestmodel))

# convergence plot
pdf(paste(opt$out,"/Convergence_ELBO",opt$suffix,".pdf",sep=""),h=4,w=4*2)
par(mfrow=c(1,1),las=1,family="Times")
plot( diff(MOFAobject@TrainStats$elbo),type="l",xlab="Iteration",ylab=expression(Delta~"ELBO"),main="",log="y")
abline(h=0,lty=2)
dev.off()

# analyses
# scatter plot 
MOFAfactors <- MOFAobject@Expectations$Z
kml = lapply(2:4, function(i) kmeans(MOFAfactors[,c(2,3)],centers = i) )
#km3 = kmeans(MOFAfactors[,c(2,3)],centers = 3)
#km2 = kmeans(MOFAfactors[,c(2,3)],centers = 2)

require(cluster)
lsil = vector("list",(4-1))
for(i in 2:4){
  sil = silhouette(kml[[i-1]]$cluster,dist(MOFAfactors[,c(2,3)],method = "euclidean"))
  sizes = table(kml[[i-1]]$cluster)
  plot( sil ,col=rep( rainbow(i),rep=sizes) ,main=paste("K=",i))
  lsil[[i-1]]=sil
}
msil = sapply(1:(4-1), function(i) mean( lsil[[i]][,3] ) )
kmb = kml[[which.max(msil)]]
K = length( kmb$size )

svg( paste(opt$out,"/ScatterPlot",opt$suffix,".svg",sep=""),h=4,w=4)
par(mfrow=c(1,1),family="Times",las=1)
plot(MOFAfactors[,2:3], xlab = "Latent Factor 1", ylab="Latent Factor 2", col =rainbow(length(kmb$size))[as.factor(kmb$cluster)], pch=16 )
legend("topleft",legend = paste("Cluster ",1:length(kmb$size),sep="") , col = rainbow(length(kmb$size)), pch=16 )
dev.off()

pdf( paste(opt$out,"/ScatterPlot",opt$suffix,"_text.pdf",sep=""),h=12,w=12)
par(mfrow=c(1,1),family="Times",las=1)
plot(-100,-100,xlim=range(MOFAfactors[,2]),ylim=range(MOFAfactors[,3]) , xlab = "Latent Factor 1", ylab="Latent Factor 2")
text(MOFAfactors[,2:3], labels =  rownames(MOFAfactors), col =rainbow(length(kmb$size))[as.factor(kmb$cluster)], pch=16 )
legend("topleft",legend = paste("Cluster ",1:length(kmb$size),sep="") , col = rainbow(length(kmb$size)), pch=16 )
dev.off()

# correlations
require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

viewsM = names( MOFAobject@TrainData )
corlist = vector("list",length(viewsM))
LF = 1:3
for(i in 1:length(viewsM)){
    print(i)
    corlist[[i]] = vector("list",length(LF))
    for(j in LF){
      print(j)
      cor.Mut = apply( MOFAobject@TrainData[[i]], 1, function(x) unlist(cor.test(x , MOFAfactors[,j+1])[c("estimate","p.value")]) )
      cor.Mut = rbind(cor.Mut, p.value.adj = p.adjust(cor.Mut[2,],method = "BH"))
      if( length(which( (cor.Mut[3,]<0.05) ))>1 ) cor.Mut = cor.Mut[,which( (cor.Mut[3,]<0.05) ) ]
      else cor.Mut = matrix( cor.Mut[,which( (cor.Mut[3,]<0.05) ) ] ,3,1,dimnames = list( rownames(cor.Mut),rownames(cor.Mut)[which( (cor.Mut[3,]<0.05) )]  ))
      print(cor.Mut)
      cor.Mut = cor.Mut[,order(abs(cor.Mut[1,]),decreasing = T )]
      if(length(grep("Meth", viewsM[i]))>0 ){
        cor.Mut = rbind(cor.Mut, t(as.matrix(ann850k[ann850k$Name%in%colnames(cor.Mut),c(22,24,31)])) )
      }
      write.table(t(cor.Mut), file = paste(opt$out,"/cor_MOFA",opt$suffix,"_LF",j,"_",viewsM[i],".txt",sep=""),quote = F)
      corlist[[i]][[j]] = cor.Mut
    }
}

# correlate other variables
if(!is.null(opt$group_file)){
  meta = read.table(opt$group_file,h=T,sep="\t")

  idsurv = grep("surv", colnames(meta) ,ignore.case = T)
  if( length(idsurv)>0){
    surv = meta[,idsurv]
    meta = meta[,-idsurv]
    issurv = TRUE
  }else{
    issurv = FALSE
  }
  
  MOFAfactors.nona = MOFAfactors[ rowSums( is.na( meta ) )==0 ,]
  
  lml = lapply(1:3, function(i) lm( as.formula(paste("MOFAfactors.nona[,",i+1,"]~", paste(colnames(meta) ,collapse = "+"),"+ Sex:Type + Sex:Smoking_Hx + Age:Smoking_Hx + Smoking_Hx:Type") ),data = na.omit(meta) ) )
  suml = lapply(lml,summary)
  anl = lapply(lml,anova)
  
  require(MASS)
  #lml.dropped = vector("list",3)
  #for(i in 1:3){
  #  lmtmp = lml[[i]]
  #  metatmp = meta
  #  while(1){
  #    dttmp = dropterm(lmtmp)
  #    if(max(dttmp$AIC,na.rm = T)< dttmp$AIC[1]  ) ttmp = rownames(dttmp)[which.min(dttmp$AIC)]
  #    else break
      #if(ttmp=="<none>") break
  #    metatmp = metatmp[,!(colnames(metatmp)==ttmp)]
  #    lmtmp = lm( as.formula(paste("MOFAfactors[,",i+1,"]~", paste(colnames(metatmp),collapse = "+")) ),data = metatmp ) 
      #lmtmp2 = update(lmtmp, paste(".~. - ",ttmp))
      #if( anova(lmtmp,lmtmp2)$`Pr(>F)`[2] >0.05/3) lmtmp = lmtmp2
  #    print(dttmp)
  #    #print(lmtmp)
  #  }
  #  lml.dropped[[i]] = lmtmp
  #}
  #anl.dropped = lapply(lml.dropped,anova,)
  lml.dropped = lapply(1:3, function(i) stepAIC(lml[[i]], direction = "both",test="F") )
  lml.final = lapply(1:3, function(i) lm( as.formula(paste("MOFAfactors[,",i+1,"]~", as.character(lml.dropped[[i]]$call$formula)[3] ) ) ,data = meta) )
  for(i in 1:3) write.table( anl[[i]],file = paste(opt$out,"/Anova",opt$suffix,"_LF",i,".txt",sep=""),quote = F) 
  anl.final = lapply(lml.final,anova)
  
  lsignif = lapply(1:3, function(i) rownames(anl.final[[i]])[which( anl.final[[i]]$`Pr(>F)`<0.05/3 )] )
  #lsignifall = unique( unlist(strsplit(unlist( lsignif ),":")) )
  lsignifall = unique( unlist(lsignif) )
  
  meta2 = cbind(meta, "Type:Sex" = interaction(meta$Type,meta$Sex), "Sex:Smoking_Hx" = interaction(meta$Sex,meta$Smoking_Hx) , "Age:Smoking_Hx" = interaction(meta$Age,meta$Smoking_Hx) )
  
  svg(paste(opt$out,"/Pie",opt$suffix,"_signif.svg",sep=""),h=3*length(lsignifall),w=3*K)
  par(mfrow=c(length(lsignifall),K),family="Times",las=1)
  for(j in 1:length(lsignifall)){
    for(i in 1:K) pie( table(factor(meta2[,lsignifall[j]])[kmb$cluster==i] ) ,radius = sqrt(sum(kmb$cluster==i)/70), init.angle = 90,col= rainbow(length(levels(factor(meta2[,lsignifall[j]])))) )
  }
  dev.off()
  
  fisher.test(table(meta2[,lsignifall[j]],kmb$cluster ))
  
  survS <- Surv( as.numeric(surv[,1]), as.numeric(factor(surv[,2],levels=c("alive","dead") ))-1)
  fitS       <- survfit(survS ~ MOFAfactors[,2] + MOFAfactors[,3] )
  fitS.clust <- survfit(survS ~ factor(kmb$cluster)  )
  fitS.Type  <- survfit(survS ~ factor(meta$Type)  )
  coxS       <- coxph(survS ~ MOFAfactors[,2] + MOFAfactors[,3]) 
  coxS.clust <- coxph(survS ~ factor(kmb$cluster)  )
  coxS.Type   <- coxph(survS ~ factor(meta$Type)  )
  summary(coxS)
  summary(coxS.clust)
  summary(coxS.Type)
  extractAIC(coxS) # better
  extractAIC(coxS.clust)
  extractAIC(coxS.Type)
  #summary(fitCarc.clust2)$table
  
  svg(paste(opt$out,"/Surv",opt$suffix,".svg",sep=""),h=4,w=4)
  par(mfrow=c(1,1),las=1,family="Times")
  plot(fitS.clust,col=1:4)
  legend("bottomright",legend = levels(as.factor(kmb$cluster)),col=1:3 , lty=1)
  legend("bottomleft",legend = paste("logrank P =", format(summary(coxS.clust)$logtest[3],digits = 3,scientific = T))  )
  #plot(fitCarc.Type,col=2:1)
  #legend("bottomright",legend = levels(as.factor(TypeCarc)),col=2:1 , lty=1)
  #legend("bottomleft",legend = paste("logrank P =", format(summary(coxCarc.Type)$logtest[3],digits = 3,scientific = T))  )
  dev.off()
  
  svg(paste(opt$out,"/Surv",opt$suffix,"_Type.svg",sep=""),h=4,w=4)
  par(mfrow=c(1,1),las=1,family="Times")
  plot(fitS.Type,col=1:4)
  legend("bottomright",legend = levels(as.factor(kmb$cluster)),col=1:3 , lty=1)
  legend("bottomleft",legend = paste("logrank P =", format(summary(coxS.Type)$logtest[3],digits = 3,scientific = T))  )
  #plot(fitCarc.Type,col=2:1)
  #legend("bottomright",legend = levels(as.factor(TypeCarc)),col=2:1 , lty=1)
  #legend("bottomleft",legend = paste("logrank P =", format(summary(coxCarc.Type)$logtest[3],digits = 3,scientific = T))  )
  dev.off()
  
  if(sum(colnames(meta)=="Type")>0){
    svg( paste(opt$out,"/ScatterPlot",opt$suffix,"_Type.svg",sep=""),h=4,w=4)
    par(mfrow=c(1,1),family="Times",las=1)
    plot(MOFAfactors[,2:3], xlab = "Latent Factor 1", ylab="Latent Factor 2", col = prettycolors[c(2,3,1,6)][addNA(as.factor(meta$Type))], pch=16 )
    legend("topleft",legend = levels(addNA(as.factor(meta$Type))) , col = prettycolors[c(2,3,1,6)], pch=16 )
    dev.off()
  }
  
  svg(paste(opt$out,"Ki67_AC-TC-LCNEC.svg",sep=""),h=4,w=4*2)
  par(mfrow=c(1,2),las=1,family="Times")
  boxplot( MOFAobject@TrainData$RNA[which( rownames( MOFAobject@TrainData$RNA )=="MKI67" ),] ~ meta$Type  )
  boxplot( MOFAobject@TrainData$RNA[which( rownames( MOFAobject@TrainData$RNA )=="MKI67" ),] ~ kmb$cluster  )
  dev.off()
  
}

save(MOFAobject,MOFAfactors, file = paste(opt$out,"/integration_MOFA",opt$suffix,".RData",sep="") )
