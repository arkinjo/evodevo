#denv <- c() 
#Plali <- c()
#PC1ali <- c()
#norm.dp <- c()
#modeltype <- c()
#sv1 <- c()
#psv1 <- c()

#envmode <- 2 ### 0: Ancestral, 1: Novel, 2: Diff
pert <- "pG" ### "pG": Genetic/Mutational perturbation, "pe": environmental perturbation

#Reminder that R starts counting at 1 instead of 0.
pertname <- switch(pert, "pG" = "Mutation", "pe" = "Environmental Change") 

df.ali100 <- data.frame(id=c(1:20))
df.sv1100 <- data.frame(id=c(1:20))
df.psv1100 <- data.frame(id=c(1:20))
df.Fnorm <- data.frame(id=c(1:20))
#Challenge: Try to plot all models in one plot.

for (layers in c("EFGHJP","CEFGHJP")){
  for (envmode in c(0,1)){ 
    modelname <- switch(layers, "EFGHJP"="Full", "_FGHJP"="NoCue", "E_G__P"="NoHier", "EFGH__"="NoDev", "__G___"="Null", "CEFGHJP"="NoTrain"  )
    modelenv <- paste(modelname,envmode)
    
    envmodename <- switch(envmode+1, "Ancestral", "Novel", "Novel-Ancestral") 
    
    ali_in <- read.table(paste(layers,"_",envmode,"_",pert,"_ali.dat",sep=""),header=FALSE)
    df.ali100[modelenv] <- ali_in$V2[which(ali_in$V1==100)] #Alignment between first singular vector with environmental change
    
    sv1_in <- read.table(paste(layers,"_",envmode,"_",pert,"_sval1.dat",sep=""),header=FALSE)
    df.sv1100[modelenv] <- sv1_in$V2[which(sv1_in$V1==100)] #First singular value
    df.psv1100[modelenv] <- sv1_in$V3[which(sv1_in$V1==100)] #Proportion first singular value
    
    Fnorm_in <- read.table(paste(layers,"_",envmode,"_",pert,"_denv22.dat",sep=""),header=FALSE)
    df.Fnorm[modelenv] <- Fnorm_in$V2[which(sv1_in$V1==100)] #Frobenius norm.
    
    
  }
}

## Colors now indicate ancestral or novel environment
cols <- c("magenta","cyan")
colvec <- rep.int(cols,times=4)
boxatvec <- c(1,2,4,5)
#boxatvec <- boxatvec[-3*seq(1,4)]
axisatvec <- c(1.5,4.5)

#tiff("derr_AR.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
tiff(sprintf("ali_%s-train.tif",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.ali100[2:ncol(df.ali100)], ylab = "Alignment", col=cols,ylim=c(0,1), at=boxatvec, xaxt = "n",cex.lab=1.5)
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("topright", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
#legend("bottomleft",legend = c("Full","NoHier","NoCue","NoDev"), col=c("orange","navy","limegreen","darkorchid"), lty=c(rep(1,4)))
dev.off()

tiff(sprintf("sv1_%s-train.tif",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.sv1100[2:ncol(df.sv1100)], xlab="Model", ylab = "SV1", col=colvec, at= boxatvec, xaxt="n")
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("topleft", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev"), col=c("orange","navy","limegreen","darkorchid"), lty=c(rep(1,4)))
dev.off()

tiff(sprintf("psv1_%s-train.tif",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.psv1100[2:ncol(df.psv1100)], xlab="Model", ylab ="% 1st singular value", col=colvec, ylim=c(0,1), at = boxatvec, xaxt = "n",cex.lab=1.5)
axis(side = 1, at = axisatvec, labels=c("Full","NoTrain"))
legend("topleft", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev"), col=c("orange","navy","limegreen","darkorchid"), lty=c(rep(1,4)))
dev.off()

tiff("Fnorm-train.tif",width=2250,height=2250,units="px",pointsize=12,res=300)
par(mar=c(5,4,4,4)+0.1)
#boxplot(df.Fnorm[2:5], main=paste(envmodename, pertname, sep=","), xlab="Model", ylab="Square Frobenius Norm", col=c("orange","limegreen","navy","darkorchid"), ylim=c(0,2.5))
boxplot(df.Fnorm[c(2,3,4,5)], ylab="Total Cross Covariance (Pheno-Cue)", col=cols, at = c(1,2,4,5), xaxt = "n", xlim=c(0.0,12.0),cex.lab=1.5) #May need to adjust range manually
axis(side = 1, at = c(1.5,4.5), labels=c("Full (Pheno-Cue)","NoTrain (Pheno-Geno)"))
par(new=TRUE)
boxplot(df.Fnorm[c(6,7,8,9)], col=cols,  at = c(7,8,10,11), xaxt = "n", xlim=c(0,12),cex.lab=1.5) #May need to adjust range manually

#legend("topleft", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
dev.off()