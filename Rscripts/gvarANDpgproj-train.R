
###Readme
#Run this after running gvar-train.R and PGtraj-train.R
#Set working directory to proj file directory before running this file

df.approj_mingvar <- data.frame(epoch=seq(1,nepoch))
df.npproj_mingvar <- data.frame(epoch=seq(1,nepoch))
df.pproj_mingvar <- data.frame(epoch=seq(1,nepoch))

df.vars_mingvar <- data.frame(epoch=seq(1,nepoch))

df.gproj_mingvar <- data.frame(epoch=seq(1,nepoch))
df.apvar_mingvar <- data.frame(epoch=seq(1,nepoch))
df.npvar_mingvar <- data.frame(epoch=seq(1,nepoch))

df.gvar_mingvar <- data.frame(epoch=seq(1,nepoch))

for(i in seq(2,ncol(df.tminvar))){
  aphat <- c()
  nphat <- c()
  ghat <- c()
  var_ghat <- c()
  var_aphat <- c()
  var_nphat <- c()
  
  model <- colnames(df.tminvar)[i]
  modelgens <- df.tminvar[model]
  
  layers <- switch(model,"Full"="EFGHJP","NoCue"="_FGHJP","NoHier"="E_G__P","NoDev"="EFGH__", "NoTrain"="CEFGHJP")
  for(j in seq(1,nrow(df.tminvar))){
    #modelXepoch <- sprintf("%s_%02d",model,j)
    gen_minvar <- modelgens[j,]
    filename <- sprintf("%s_run100_%02d_%03d.dat",layers,j,gen_minvar)
    pgstats <- read.delim(filename,skip = 1000)
    
    aphat <- append(aphat,pgstats[1,4])
    var_aphat <- append(var_aphat,pgstats[1,6]^2)
    nphat <- append(nphat,pgstats[2,4])
    var_nphat <- append(var_nphat,pgstats[2,6]^2)
    
    ghat <- append(ghat,pgstats[2,3])
    var_ghat <- append(var_ghat,pgstats[2,5]^2)
  }
  df.approj_mingvar[model] <- aphat
  df.npproj_mingvar[model] <- nphat
  df.pproj_mingvar[sprintf("%s_a",model)] <- aphat
  df.pproj_mingvar[sprintf("%s_n",model)] <- nphat
  
  df.vars_mingvar[sprintf("%s_a",model)] <- var_aphat
  df.vars_mingvar[sprintf("%s_n",model)] <- var_nphat
  
  df.gproj_mingvar[model] <- ghat
  df.gvar_mingvar[model] <- var_ghat
  
}

# boxatvec <- seq(1,12)
# boxatvec <- boxatvec[-3*seq(1,4)]
# axisatvec <- 3*seq(1,4)-1.5

boxatvec <- c(1,2,4,5)
axisatvec <- c(1.5,4.5)

tiff("../plots/approjhat-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.approj_mingvar[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="Projected Novel Phenotype",ylim=c(0,1))
dev.off()

tiff("../plots/npprojhat-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.npproj_mingvar[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="Projected Novel Phenotype",ylim=c(0,1))
dev.off()

tiff("../plots/pprojhat-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.pproj_mingvar[seq(2,ncol(df.pproj_mingvar))],col=c("darkorchid","cyan"),at=boxatvec,xlab="Model",ylab="Projected Phenotype",main="Genetic bottleneck",ylim=c(0,1),xaxt="n")
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("bottomright",legend=c("Ancestral","Novel"),col=c("darkorchid","cyan"),lty=1,title="Environment")
dev.off()

tiff("../plots/varpprojhat-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.vars_mingvar[seq(2,ncol(df.vars_mingvar))],col=c("darkorchid","cyan"),at=boxatvec,xlab="Model",ylab="Variance in projected phenotype",xaxt="n")
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("topright",legend=c("Ancestral","Novel"),col=c("darkorchid","cyan"),lty=1,title="Environment")
dev.off()


tiff("../plots/gprojhat-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.gproj_mingvar[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="Projected Genotype",ylim=c(0,1))
dev.off()

#tiff("../plots/gvarhat.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
#boxplot(df.gvar_mingvar[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="Decrease in error",main="Change in regulation")
#dev.off()
