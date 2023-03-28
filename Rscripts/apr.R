setwd("../all")

#maxgen <- 200 #Epoch length
nepoch <- 10
popsize = 1000
gen = 200

#dG <- data.frame(gen=c(1:(maxgen-1)))
#npvar <- data.frame(gen=c(1:maxgen))
#apvar <- data.frame(gen=c(1:maxgen))

df.ipr <- data.frame(id=c(1:popsize))
df.mipr <- data.frame(id = c(1:(popsize*nepoch)))
df.vars <- data.frame(id=c(1:10)) #Put everything in same plot


#df.sdap <- data.frame(id=c(1:10))
#df.sdnp <- data.frame(id=c(1:10))

models <- c()
colvec.s <- c()
colvec.l <- c()
#modelipr <- c()

for(layers in c("EFGHJP","E_G__P","_FGHJP", "EFGH__")){ #Loop over models
  
  modelipr <- c() #Empty for each new model
  sdap <- c() #Standard deviation in ancestral phenotype
  sdnp <- c() #Standard deviation in novel phenotype
  
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","__G___"="Null")
  models <- append(models,modelname)
  modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid")
  colvec.s <- append(colvec.s,modelcol)
  
  for(epoch in c(1:nepoch)){
    str.epoch <- sprintf("%02d",epoch)
    t0gen <- sprintf("%03d",gen)
    filename <- paste(layers,"_run100_",str.epoch,"_",t0gen,".dat",sep="")
    pgstats <- read.delim(filename,nrows = 1000)
    #View(pgstats)
    ipr <- pgstats$Pheno1-pgstats$Pheno0
    df.ipr[paste(modelname,str.epoch)] <- ipr
    modelipr <- append(modelipr, ipr)
    
    avesd <- read.delim(filename,skip=1000)
    sdap <- append(sdap,avesd[1,6])
    sdnp <- append(sdnp,avesd[2,6])
    
    colvec.l <- append(colvec.l, modelcol)
  }
  
  df.mipr[modelname] <- modelipr
  df.vars[paste(modelname,"ap",sep="_")] <- sdap^2
  df.vars[paste(modelname,"np",sep="_")] <- sdnp^2
  
}

atvec <- c(1:12)
atvec <- atvec[-c(1:4)*3]


tiff(filename=sprintf("../plots/apr_%d.tif",gen),width=2250,height=2250,units="px",pointsize=12,res=300)
bp.ipr <- boxplot(df.ipr[2:length(df.ipr)],col=colvec.l,las=2,ylab="Adaptive plastic response",main=sprintf("Generation %d",gen),ylim=c(-1,1))
#legend("topright",legend=c("Full","NoHier","NoDev"),col=c("orange","limegreen","darkorchid"),lty=c(1,1,1))
dev.off()

tiff(filename=sprintf("../plots/mapr_%d.tif",gen),width=2250,height=2250,units="px",pointsize=12,res=300)
bp.mipr <- boxplot(df.mipr[2:length(df.mipr)],col=colvec.s,las=1, xlab="Model", ylab="Adaptive plastic response", main=sprintf("Generation %d",gen))
#legend("topright",legend=c("Full","NoHier","NoDev"),col=c("orange","limegreen","darkorchid"),lty=c(1,1,1))
dev.off()

tiff(filename=sprintf("../plots/vars_%d.tif",gen),width=2250,height=2250,units="px",pointsize=12,res=300)
bp.vars <- boxplot(df.vars[2:length(df.vars)],col=c("magenta","cyan"),ylab=sprintf("Variance in projected phenotype (Gen %d)",gen), at = atvec, xaxt="n",cex.lab=1.5,cex.main=2.0)
axis(side = 1, at = c(1.5,4.5,7.5,10.5), labels = models )
legend("topright",legend=c("Ancestral", "Novel"), col=c("magenta","cyan"), lty=1, title="Environment")
dev.off()

#png(filename=paste("../plots/sdap",gen,".png"))
#bp.sdap <- boxplot(df.sdap[2:4],col=c("orange","limegreen","darkorchid"),las=1, ylab="Standard deviation in phenotype", main=paste("Generation",gen,sep=" "))
#legend("topright",legend=c("Full","NoHier","NoDev"),col=c("orange","limegreen","darkorchid"),lty=c(1,1,1))
#dev.off()

#png(filename=paste("../plots/sdnp",gen,".png"))
#bp.sdnp <- boxplot(df.sdnp[2:4],col=c("orange","limegreen","darkorchid"),las=1, ylab="Standard deviation in phenotype", main=paste("Generation",gen,sep=" "))
#legend("topright",legend=c("Full","NoHier","NoDev"),col=c("orange","limegreen","darkorchid"),lty=c(1,1,1))
#dev.off()

# Challenge: Put ancestral and novel of each model side by side, all boxplots for cryptic mutations in one figure. 
# Will need to adjust boxplot axis parameters.
