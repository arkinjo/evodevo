
nepoch <- 10
ngen <- 200
fpt.T <- 0.20 #Arbitrary number


ffpt.s <- function(v,Omega) { #First passage time
  return(min(which(v<Omega)))
}

df.fpt <- data.frame(epoch=1:nepoch)
df.err <- data.frame(gen=1:ngen)
df.err_inf <- data.frame(epoch=seq(1,nepoch))
#df.diff.err <- data.frame(gen=c(1:(ngen-1)))
df.pdiv <- data.frame(gen=1:ngen)


colvec.s <- c() #Boxplot
colvec.l <- c() #Trajectory
modelvec <- c() #Legend

for(layers in c("EFGHJP","CEFGHJP")){ #Loop over models, order determines what goes over/under in plot
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","__G___"="Null","CEFGHJP"="NoTrain")
  traj <- read.table(sprintf("%s_run100.traj",layers),header=FALSE)
  modelvec <- append(modelvec, modelname)
  modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid", "__G___"="red","CEFGHJP"="blue")
  colvec.s <- append(colvec.s, modelcol)
  modelfpt <- c()
  err_infty <- c()
  for (epoch in c(1:nepoch)){
    modelXepoch <- sprintf("%s_%02d",modelname,epoch)
    errtraj <- traj$V5[which(traj$V1==epoch)]
    pdivtraj <- traj$V12[which(traj$V1==epoch)]
    df.err[modelXepoch] <- errtraj
    #df.diff.err[modelXepoch] <- diff(errtraj)
    df.pdiv[modelXepoch] <- pdivtraj
    err_infty <- append(err_infty,errtraj[ngen])
    modelfpt <- append(modelfpt,ffpt.s(errtraj,fpt.T))
    colvec.l <- append(colvec.l,modelcol)
    
  }
  
  df.err_inf[modelname] <- err_infty
  df.fpt[modelname] <- modelfpt
}


bxpindex <- c(1,2)
mat.err <- as.matrix(df.err)
#mat.diff.err <- as.matrix(df.diff.err)
mat.pdiv <- as.matrix(df.pdiv)


#png("pdivtraj.png")
tiff("pdivtraj-train.tif",width=2250,height=2250,units="px",pointsize=12,res=300)
matplot(mat.pdiv[,2:ncol(mat.pdiv)],col=colvec.l,type="l",lty=1,xlab="Generation",ylab="Phenotypic Variance")
legend("topright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
dev.off()



tiff("errtraj-train.tif",width=2250,height=2250,units="px",pointsize=12,res=300)
matplot(mat.err[,2:ncol(mat.err)],col=colvec.l,type="l",lty=1,xlab="Generation",ylab="Absolute error")
legend("topright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
dev.off()

tiff("err_infty-train.tif",width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.err_inf[bxpindex+1],col=colvec.s[bxpindex],type="l",lty=1,xlab="Model",ylab="Final error")
#legend("topright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
dev.off()

#png("differrtraj.png")
#matplot(mat.diff.err[,2:ncol(mat.diff.err)],col=colvec.l,type="l",lty=1,xlab="Generation",ylab="Change in error")
#abline(h=0,col="black")
#legend("bottomright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
#dev.off()


#png("errfpt.png")
#boxplot(df.fpt[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="First passage time",main="Absolute error")
#legend("topleft",legend=modelvec[bxpindex], col=colvec.s[bxpindex],lty=1)
#dev.off()
