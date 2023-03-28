
#layerlist <- c("EFGHJP","E_G__P","_FGHJP")

#Make sure to set working directory to source file location prior to 
setwd("../all")

maxgen <- 200 #Epoch length
nepoch <- 10
#fpt.Threshold <- 0.9 #Arbitrary number between 0 and 1.
#str.fpt <- sprintf("%.2f",fpt.Threshold)

df.G <- data.frame(gen=c(1:maxgen))
df.dG <- data.frame(gen=c(1:(maxgen-1)))
npvar <- data.frame(gen=c(1:maxgen))
apvar <- data.frame(gen=c(1:maxgen))

df.pGvar <- data.frame(gen=c(1:maxgen))

#ffpt.s <- function(v,Omega) { #First passage time
#  return(min(which(v<Omega)))
#}

#ffpt.l <- function(v,Omega) { #First passage time
#  return(min(which(v>Omega)))
#}

#No need first passage times
#df.fpt.ag <- data.frame(epoch = c(1:nepoch))
#df.fpt.ap <- data.frame(epoch = c(1:nepoch)) #First passage time for this is not sane
#df.fpt.ng <- data.frame(epoch = c(1:nepoch))
#df.fpt.np <- data.frame(epoch = c(1:nepoch))

df.p_0 <- data.frame(epoch=seq(1,nepoch)) #Projected phenotypes before evolution
df.p_inf <- data.frame(epoch=seq(1,nepoch)) #Projected phenotypes after evolution


#models <- c("Full","NoCue","NoHier","NoDev")
modelvec <- c()
colvec.s <- c()
colvec.l <- c()

for(layers in c("EFGHJP","CEFGHJP")){ #Loop over models, order determines what goes over/under in plot!
  
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","__G___"="Null","CEFGHJP"="NoTrain")
  modelvec <- append(modelvec, modelname)
  modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid", "__G___"="red","CEFGHJP"="blue")
  colvec.s <- append(colvec.s, modelcol)
  
  ap_0 <- c()
  np_0 <- c()
  ap_inf <- c()
  np_inf <- c()
  
  #print(modelname)
  #fpt.ag <- c()
  #fpt.ap <- c()
  #fpt.np <- c()
  #fpt.ng <- c()
  
  #all.nghat <- c() 
  #all.npvar <- c()
  #all.apvar <- c()
  
  for(epoch in c(1:nepoch)) {
    
    
    str.epoch <- sprintf("%02d",epoch)
    colvec.l <- append(colvec.l,modelcol)
    
    aghat <- c()
    aphat <- c()
    nghat <- c()
    nphat <- c()
    
    gvar <- c()
    
    #fpt.aghat <- 0
    #fpt.aphat <- 0
    #fpt.nphat <- 0
    #fpt.nghat <- 0
    
    sdap <- c()
    sdnp <- c()
    sdag <- c()
    sdng <- c()
    
    #pdffilename <- paste("nodev.pdf",sep="")
    #pdf(pdffilename)
    for (gen in c(1:maxgen)){
      t0gen <- sprintf("%03d",gen)
      filename <- paste(layers,"_run100_",str.epoch,"_",t0gen,".dat",sep="")
      
      pgstats <- read.delim(filename,skip = 1000)
      #print(filename)
      
      aghat <- append(aghat,pgstats[1,3])
      aphat <- append(aphat,pgstats[1,4])
      sdag <- append(sdag,pgstats[1,5])
      sdap <- append(sdap,pgstats[1,6])
      nghat <- append(nghat,pgstats[2,3])
      nphat <- append(nphat,pgstats[2,4])
      sdng <- append(sdng,pgstats[2,5])
      sdnp <- append(sdnp,pgstats[2,6])
    }
    #dev.off()
    
    #Summarize trajectory: error plots for initial, reference and final generation.
    #pdf(paste("../plots/",modelname,"_",str.epoch,".pdf",sep=""))
    
    
    tiff(sprintf("../plots/%s_%02d.tif",modelname,epoch),width=2250,height=2250,units="px", pointsize=12, res=300)
    plot(aghat,aphat,xlim=c(0.0,1.0),ylim=c(0.0,1.0),xlab="Genotype",ylab="Phenotype",main=modelname,col="darkorchid",type="l",cex.lab=1.5,cex.main=2.0)
    points(aghat,aphat,col="darkorchid")
    arrows(aghat, aphat-sdap, aghat, aphat+sdap, length=0.05, angle=90, code=3, col="darkorchid") 
    arrows(aghat-sdag, aphat, aghat+sdag, aphat, length=0.05, angle=90, code=3, col="darkorchid") 
    lines(nghat,nphat,col="cyan")
    points(nghat,nphat,col="cyan")
    arrows(nghat, nphat-sdnp, nghat, nphat+sdnp, length=0.05, angle=90, code=3, col="cyan") 
    arrows(nghat-sdng, nphat, nghat+sdng, nphat, length=0.05, angle=90, code=3, col="cyan") 
    legend("topleft",legend=c("Novel","Ancestral"),col=c("cyan","darkorchid"),title="Environment",lty=c(1,1))
    dev.off()
    
    df.G[paste(modelname,str.epoch,sep="")] <- nghat
    df.dG[paste(modelname,str.epoch,sep="")] <- diff(nghat)
    apvar[paste(modelname,str.epoch,sep="")] <- sdap^2
    npvar[paste(modelname,str.epoch,sep="")] <- sdnp^2
    
    df.pGvar[sprintf("%s%02d",modelname,epoch)] <- sdng^2
    
    ap_0 <- append(ap_0,aphat[1])
    np_0 <- append(np_0,nphat[1])
    ap_inf <- append(ap_inf,aphat[maxgen])
    np_inf <- append(np_inf,nphat[maxgen])
    
    #ap_inf <- append(ap_inf)
    
    
    #fpt.aghat <- ffpt.l(aghat,fpt.Threshold)
    #fpt.aphat <- ffpt.l(aphat,fpt.Threshold)
    #fpt.nghat <- ffpt.l(nghat,fpt.Threshold)
    #fpt.nphat <- ffpt.l(nphat,fpt.Threshold)
    
    #fpt.ag <- append(fpt.ag, fpt.aghat)
    #fpt.ap <- append(fpt.ap, fpt.aphat)
    #fpt.ng <- append(fpt.ng, fpt.nghat)
    #fpt.np <- append(fpt.np, fpt.nphat)
    print(paste("Plotting ../plots/",modelname,"_",str.epoch,".tif",sep=""))
  }
  
  df.p_0[sprintf("%s_a",modelname)] <- ap_0
  df.p_0[sprintf("%s_n",modelname)] <- np_0
  df.p_inf[sprintf("%s_a",modelname)] <- ap_inf
  df.p_inf[sprintf("%s_n",modelname)] <- np_inf
  
  
  
  #df.fpt.ag[modelname] <- fpt.ag
  #df.fpt.ap[modelname] <- fpt.ap
  #df.fpt.ng[modelname] <- fpt.ng
  #df.fpt.np[modelname] <- fpt.np
  
}

#Change working directory to NoDev here! #No need anymore
#setwd("../nodev")

mat.G <- as.matrix(df.G)
#dG <- diff(df.G) #Bad line
mat.dG <- as.matrix(df.dG)

mat.apvar <- as.matrix(apvar)
mat.npvar <- as.matrix(npvar)

mat.pGvar <- as.matrix(df.pGvar)

#lindex <- c(3,2,1,4)
lindex <-
boxatvec <- c(1,2,4,5)
axisatvec <- c(1.5,4.5)



tiff("../plots/G_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
matplot(y=mat.G[,c(2:ncol(mat.G))],ylab="Projected genotype",xlab="Generation",type="l",lty=1,col=colvec.l)
legend("topright",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)



#pdf("../plots/dG.pdf")
tiff("../plots/dG_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
matplot(y=mat.dG[,c(2:ncol(mat.dG))],ylab="Change in projected genotype",xlab="Generation",type="l",lty=1,col=colvec.l)
legend("topright",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)

#legend("topright",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))

dev.off()

#pdf("../plots/dG-log.pdf")
tiff("../plots/dG-log_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)

matplot(y=mat.dG[,c(2:ncol(mat.dG))],ylab="Change in projected genotype",xlab="Generation",type="l",lty=1,col=colvec.l,log="y")
legend("topright",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)

#legend("topright",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))
dev.off()

#pdf("../plots/npvar.pdf")
tiff("../plots/npvar_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)

matplot(y=mat.npvar[,c(2:ncol(mat.npvar))],ylab="Variance in novel phenotype",xlab="Generation",type="l",lty=1,col=colvec.l)
legend("topright",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))
dev.off()

tiff("../plots/pGvar_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
matplot(y=mat.pGvar[,c(2:ncol(mat.pGvar))],ylab="Variance in projected genotype",xlab="Generation",type="l",lty=1,col=colvec.l)
legend("topright",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))
dev.off()



tiff("../plots/npvar-log_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
matplot(y=mat.npvar[,c(2:ncol(mat.npvar))],ylab="Variance in novel phenotype",xlab="Generation",type="l",lty=1,col=colvec.l,log="y")
legend("topright",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))
dev.off()

tiff("../plots/apvar_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)

matplot(y=mat.apvar[,c(2:ncol(mat.apvar))],ylab="Variance in ancestral phenotype",xlab="Generation",type="l",lty=1,col=colvec.l)
legend("topleft",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)

#legend("topright",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))
dev.off()

tiff("../plots/apvar-log_train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)

#png("../plots/apvar-log.png")
matplot(y=mat.apvar[,c(2:ncol(mat.apvar))],ylab="Variance in ancestral phenotype",xlab="Generation",type="l",lty=1,col=colvec.l,log="y")
legend("topleft",legend = modelvec[lindex], col=colvec.s[lindex], lty=1)
#legend("topleft",legend = c("Full","NoHier","NoCue","NoDev","Null"), col=c("orange","navy","limegreen","darkorchid","red"), lty=c(rep(1,4)))
dev.off()

tiff("../plots/p_0-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.p_0[,c(2:ncol(df.p_0))],col=c("magenta","cyan"),at=boxatvec,ylab="Projected Phenotype (Gen 1)",ylim=c(0,1),xaxt="n",cex.lab=1.5)
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("topleft",legend=c("Ancestral","Novel"),col=c("magenta","cyan"),lty=1,title="Environment")
dev.off()

tiff("../plots/p_inf-train.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.p_inf[2:ncol(df.p_0)],col=c("magenta","cyan"),at=boxatvec,ylab="Projected Phenotype (Gen 200)",ylim=c(0,1),xaxt="n",cex.lab=1.5)
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("bottomright",legend=c("Ancestral","Novel"),col=c("magenta","cyan"),lty=1,title="Environment")
dev.off()

