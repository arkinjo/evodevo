
nepoch <- 10
ngen <- 200
denv <- 70

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
epochlist <- c() 


#for(layers in c("_FGHJP","E_G__P","EFGHJP","EFGH__")){ #Loop over models, order determines what goes over/under in plot
for(layers in c("EFGHJP")){ #Loop over models, order determines what goes over/under in plot
  #modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","__G___"="Null")
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue")
  modelvec <- append(modelvec, modelname)
  #modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid", "__G___"="red")
  modelcol <- switch(layers, "EFGHJP" = "red", "_FGHJP"="navy")
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
    epochlist <- append(epochlist,sprintf("%02d",epoch+40))
  }
  df.gvar0[modelname] <- gvar0
  df.minvar[modelname] <- minvar
  df.tminvar[modelname] <- tminvar
  df.Dvar[modelname] <- Dvar
  df.Cvar[modelname] <- Cvar
}

#bxpindex <- c(3,2,1,4)
#bxpindex <- c(2,1)


#mat.agvar <- as.matrix(df.agvar)
mat.ngvar <- as.matrix(df.ngvar)
mat.diffgvar <- as.matrix(df.diffgvar)

#png(sprintf("ngvar_%d.png",denv),width=2250,height=2250,units="px", pointsize=12, res=300)
#matplot(mat.ngvar[,2:ncol(mat.ngvar)],col=colvec.l,type="l",lty=1,lwd=2,
#        xlab="Generation",ylab="Genetic Variance",cex.lab=1.5)
#legend("bottomright",legend=modelvec[bxpindex],title="Model",col=colvec.s[bxpindex],lty=1)
#dev.off()

png(sprintf("ngvar_%d.png",denv),width=2250,height=2250,units="px", pointsize=12, res=300)
#make main plot first
par(fig=c(0,1,0,1)) #Specify NDC coordinates for main plot
matplot(mat.ngvar[,2:ncol(mat.ngvar)], col=c(1:nepoch), type="l",lty=c(1:nepoch),lwd=2, 
        sub=sprintf("%s Environmental change",str.pdenv),
        xlab="Generation",ylab="Genetic Variance",cex.lab=1.5)
legend("bottomright",legend=epochlist,title="Epoch",col=c(1:nepoch),lty=c(1:nepoch))
par(fig=c(0.35,0.85,0.05,0.50),new=T)
matplot(mat.ngvar[1:10,2:ncol(mat.ngvar)],col=c(1:nepoch),
        type="l",lty=c(1:nepoch),xlab="",ylab="") #Remove axis labels
dev.off()


