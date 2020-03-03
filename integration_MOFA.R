####################################################
### Script to perform MOFA                       ###
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
  make_option(c("-r", "--robustness"), action="store_true", default=FALSE, help="Perform robustness analysis [default= %default]"),
  make_option(c("-i", "--maxiter"), type="numeric", default=10000, help="Maximum number of iterations [default= %default]"),
  make_option(c("-n", "--nconsreps"), type="numeric", default=100, help="Number of iterations for consensus clustering [default= %default]"),
  make_option(c("-m", "--nreps"), type="numeric", default=10, help="Maximum number of iterations for Latent Factor robustness analysis [default= %default]"),
  make_option(c("-d", "--dropLFthres"), type="numeric", default=0, help="Threshold proportion of variance explained to drop latent factors [default= %default]"),
  make_option(c("-t", "--tol"), type="numeric", default=0.01, help="Tolerated delta ELBO to call convergence [default= %default]"),
  make_option(c("-F", "--nLF"), type="numeric", default=10, help="Number of latent factors [default= %default]"),
  make_option(c("-p", "--sparsity"), type="logical", action="store_false", default=TRUE, help="Sparsity of latent factors [default= %default]")
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
    samptmp = sort(unique( unlist(sapply(data,colnames)) ))
    nsamp = length(samptmp)
    samptmp = sort(sample(samptmp,size = round(pItem*nsamp),replace=F) )
    data2 = lapply(data, function(x) x[,unlist(sapply(samptmp, function(xx) which(colnames(x)==xx) )) ] )
    print(sapply(data2,dim))
    V2    = lapply(data2, function(x) apply(x,1,var) )

    zeroV2 = lapply(V2, function(x) which(x==0))
    for(j in 1:length(data2)){
	print(paste("Remove",length(zeroV2[[j]]),"samples with 0 variance in view",j) )
	if( length(zeroV2[[j]])>0 ) data2[[j]] = data2[[j]][-zeroV2[[j]],]
    }
    print(sapply(data2,dim))    
    MOFAobjecttmp <- createMOFAobject(data2)  
    ModelOptions <- getDefaultModelOptions(MOFAobjecttmp)
    ModelOptions$likelihood[names(ModelOptions$likelihood)=="Mutation"] = "bernoulli"
    ModelOptions$numFactors  <- as.numeric(opt$nLF)
    ModelOptions$sparsity    <- as.logical(opt$sparsity)
    TrainOptions <- getDefaultTrainOptions()
    TrainOptions$maxiter = as.numeric(opt$maxiter)
    TrainOptions$DropFactorThreshold = as.numeric(opt$dropLFthres)
    TrainOptions$tolerance = as.numeric(opt$tol)
    DataOptions <- getDefaultDataOptions()
    
    print(ModelOptions)
    print(TrainOptions)
    print(DataOptions)
    MOFAobjecttmp <- prepareMOFA(MOFAobjecttmp, ModelOptions = ModelOptions,TrainOptions = TrainOptions,DataOptions = DataOptions)
  
    # run
    MOFAobjecttmp <- runMOFA(MOFAobjecttmp, paste(opt$out,"/MOFA",opt$suffix,"_sub",i,".hdf5",sep="") )
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
    ModelOptions <- getDefaultModelOptions(MOFAobjecttmp)
    ModelOptions$likelihood[names(ModelOptions$likelihood)=="Mutation"] = "bernoulli"
    ModelOptions$numFactors    <- as.numeric(opt$nLF)
    ModelOptions$sparsity    <- as.logical(opt$sparsity)
    TrainOptions <- getDefaultTrainOptions()
    TrainOptions$maxiter   = as.numeric(opt$maxiter)
    TrainOptions$DropFactorThreshold = as.numeric(opt$dropLFthres)
    TrainOptions$tolerance = as.numeric(opt$tol)
    DataOptions <- getDefaultDataOptions()
    print(ModelOptions)
    print(TrainOptions)
    print(DataOptions)
    MOFAobjecttmp <- prepareMOFA(MOFAobjecttmp, ModelOptions = ModelOptions,TrainOptions = TrainOptions,DataOptions = DataOptions)
    
    # run
    MOFAobjecttmp <- runMOFA(MOFAobjecttmp, paste(opt$out,"/MOFA",opt$suffix,"_run",i,".hdf5",sep="") )
    
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
