
nepoch <- 10
ngen <- 200
denv <- 10
str.pdenv <- sprintf("%d%%",denv/2)

fpt.T <- 0.20 #Arbitrary number


ffpt.s <- function(v,Omega) { #First passage time
  return(min(which(v<Omega)))
}

df.fpt <- data.frame(epoch=1:nepoch)
df.err <- data.frame(gen=1:ngen)
df.err_inf <- data.frame(epoch=1:nepoch)
#df.diff.err <- data.frame(gen=c(1:(ngen-1)))
df.derr <- data.frame(epoch=1:nepoch)
df.pdiv <- data.frame(gen=1:ngen)


colvec.s <- c() #Boxplot
colvec.l <- c() #Trajectory
modelvec <- c() #Legend

for(layers in c("_FGHJP","E_G__P","EFGHJP","EFGH__")){ #Loop over models, order determines what goes over/under in plot
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","__G___"="Null")
  traj <- read.table(sprintf("%s_run%d.traj",layers,denv),header=FALSE)
  modelvec <- append(modelvec, modelname)
  modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid", "__G___"="red")
  colvec.s <- append(colvec.s, modelcol)
  modelfpt <- c()
  err_infty <- c()
  derr <- c()
  for (epoch in c(1:nepoch)){
    modelXepoch <- sprintf("%s_%02d",modelname,epoch)
    errtraj <- traj$V5[which(traj$V1==epoch)]
    pdivtraj <- traj$V12[which(traj$V1==epoch)]
    df.err[modelXepoch] <- errtraj
    #df.diff.err[modelXepoch] <- diff(errtraj)
    df.pdiv[modelXepoch] <- pdivtraj
    
    err_infty <- append(err_infty,errtraj[ngen])
    derr <- append(derr,errtraj[1]-errtraj[ngen])
    modelfpt <- append(modelfpt,ffpt.s(errtraj,fpt.T))
    colvec.l <- append(colvec.l,modelcol)
    
  }
  
  df.err_inf[modelname] <- err_infty
  df.fpt[modelname] <- modelfpt
  df.derr[modelname] <- derr
}


bxpindex <- c(3,2,1,4)
mat.err <- as.matrix(df.err)
#mat.diff.err <- as.matrix(df.diff.err)
mat.pdiv <- as.matrix(df.pdiv)


#png("pdivtraj.png")
png("pdivtraj.png",width=2250,height=2250,units="px",pointsize=12,res=300)
matplot(mat.pdiv[,2:ncol(mat.pdiv)],col=colvec.l,type="l",lty=1,xlab="Generation",ylab="Phenotypic Variance",cex.lab=1.5)
legend("topright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
dev.off()



png(sprintf("errtraj_%d.png",denv),width=2250,height=2250,units="px",pointsize=12,res=300)
matplot(mat.err[,2:ncol(mat.err)],col=colvec.l,type="l",lty=1,
        xlab="Generation",ylab="Mismatch",sub=sprintf("(%s environmental change)",str.pdenv),cex.lab=1.5)
legend("topright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],title="Model",lty=1)
dev.off()

png("err_infty.png",width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.err_inf[bxpindex+1],col=colvec.s[bxpindex],type="l",lty=1,ylab="Mismatch (Generation 200)",cex.lab=1.5)
#legend("topright",legend=modelvec[bxpindex],col=colvec.s[bxpindex],lty=1)
dev.off()

png("derr.png",width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.derr[bxpindex+1],col=colvec.s[bxpindex],type="l",lty=1,ylab="Total decrease in mismatch",cex.lab=1.5)
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
