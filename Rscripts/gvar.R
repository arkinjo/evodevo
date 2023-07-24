
nepoch <- 10
ngen <- 200
denv <- 10

str.pdenv <- sprintf("%d%%",denv/2)
#df.agvar <- data.frame(gen=seq(1,ngen))
df.ngvar <- data.frame(gen=seq(1,ngen))
df.diffgvar <- data.frame(gen=seq(1,ngen-1))

modelvec <- c()
colvec.s <- c()
colvec.l <- c()

df.gvar0 <- data.frame(epoch=seq(1,nepoch)) #Initial value of variance
df.minvar <- data.frame(epoch=seq(1,nepoch)) #Value of minimum variance
df.tminvar <- data.frame(epoch=seq(1,nepoch)) #Generation when minimum variance happens, length of change in regulation phase
#ldf.tminvar <- data.frame(a=0) #long
#df.tminvar1 <- data.frame() #Empty data frame
df.Dvar <- data.frame(epoch=seq(1,nepoch)) #Net drop in variance
df.Cvar <- data.frame(epoch=seq(1,nepoch)) #Gain in genetic variance (accummulation of cryptic mutations)

for(layers in c("_FGHJP","E_G__P","EFGHJP","EFGH__")){ #Loop over models, order determines what goes over/under in plot
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","__G___"="Null")
  modelvec <- append(modelvec, modelname)
  modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid", "__G___"="red")
  colvec.s <- append(colvec.s, modelcol)
  
  gvar0 <- c() #Initial value of variance
  minvar <- c() #Value of minimum variance
  tminvar <- c() #Generation of minimum variance
  Dvar <- c() #Net drop in variance
  Cvar <- c()

  for (epoch in c(1:nepoch)){
    traj <- read.table(sprintf("%s_run%d_%02d.gvar",layers,denv,epoch),header=FALSE)
    modelXepoch <- sprintf("%s_%02d",modelname,epoch)
    #agvartraj <- traj$V2
    ngvartraj <- traj$V3
    #df.agvar[modelXepoch] <- agvartraj
    df.ngvar[modelXepoch] <- ngvartraj
    df.diffgvar[modelXepoch] <- diff(ngvartraj)
    gvar0 <- append(gvar0,ngvartraj[1]) #R starts counting from 1 for some reason
    minvar <- append(minvar,min(ngvartraj))
    tminvar <- append(tminvar,which.min(ngvartraj))
    Dvar <- append(Dvar,ngvartraj[1] - min(ngvartraj))
    Cvar <- append(Cvar,ngvartraj[200] - min(ngvartraj))
    #ldf.tminvar[modelXepoch] <- which.min(ngvartraj)
    #modelfpt <- append(modelfpt,ffpt.s(errtraj,fpt.T))
    colvec.l <- append(colvec.l,modelcol)
  }
  df.gvar0[modelname] <- gvar0
  df.minvar[modelname] <- minvar
  df.tminvar[modelname] <- tminvar
  df.Dvar[modelname] <- Dvar
  df.Cvar[modelname] <- Cvar
}

bxpindex <- c(3,2,1,4)
#bxpindex <- c(2,1)


#mat.agvar <- as.matrix(df.agvar)
mat.ngvar <- as.matrix(df.ngvar)
mat.diffgvar <- as.matrix(df.diffgvar)

png(sprintf("ngvar_%d.png",denv),width=2250,height=2250,units="px", pointsize=12, res=300)
matplot(mat.ngvar[,2:ncol(mat.ngvar)],col=colvec.l,type="l",lty=1,lwd=2,
        xlab="Generation",ylab="Genetic Variance",cex.lab=1.5, sub=sprintf("(%s environmental change)",str.pdenv))
legend("bottomright",legend=modelvec[bxpindex],title="Model",col=colvec.s[bxpindex],lty=1)
dev.off()

#matplot(mat.diffgvar[,2:ncol(mat.diffgvar)],col=colvec.l,type="l",lty=1,lwd=2,xlab="Generation",ylab="Genetic Variance",ylim=c(-10,10))
#legend("bottomright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
#abline(h=0)

png("minvar.tif",width=2250,height=2250,units="px",pointsize=12, res=300)
boxplot(df.minvar[bxpindex+1],col=colvec.s[bxpindex],ylab="Minimum genetic variance",cex.lab=1.5)
dev.off()

png("tminvar.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.tminvar[bxpindex+1],col=colvec.s[bxpindex],ylab="Generations to genetic bottleneck",cex.lab=1.5)
dev.off()

png("gvar0.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.gvar0[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="Initial genetic variance",cex.lab=1.5)
dev.off()

png("Dvar.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.Dvar[bxpindex+1],col=colvec.s[bxpindex],ylab="Drop in genetic variance",cex.lab=1.5)
dev.off()

png("Cvar.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.Cvar[bxpindex+1],col=colvec.s[bxpindex],ylab="Gain in genetic variance",cex.lab=1.5)
dev.off()

#df.tminvar1 <- df.tminvar[,-1]
#list.tminvar <- as.data.frame.list(df.tminvar1,row.names="0",column.names=)


#png("agvar.png")
#matplot(mat.agvar[,2:ncol(mat.agvar)],col=colvec.l,type="l",lty=1,xlab="Generation",ylab="Genetic Variance")
#legend("bottomright",legend=modelvec[bxpindex-1],col=colvec.s[bxpindex-1],lty=1)
#dev.off()

