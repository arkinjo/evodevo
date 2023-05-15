#denv <- c() 
#Plali <- c()
#PC1ali <- c()
#norm.dp <- c()
#modeltype <- c()
#sv1 <- c()
#psv1 <- c()

#envmode <- 2 ### 0: Ancestral, 1: Novel, 2: Diff
#par(cex=1.0, cex.axis=1.0, cex.lab=1.5, cex.main=2.0) #Graphical parameters #This doesn't work with boxplot for some reason

pert <- "pe" ### "pG": Genetic/Mutational perturbation, "pe": environmental perturbation

#Reminder that R starts counting at 1 instead of 0.
pert_vs <- switch(pert, "pG" = "Pheno-Geno", "pe" = "Pheno-Cue") 



df.ali100 <- data.frame(id=c(1:20))
df.sv1100 <- data.frame(id=c(1:20))
df.psv1100 <- data.frame(id=c(1:20))
df.Fnorm <- data.frame(id=c(1:20))
#Challenge: Try to plot all models in one plot.

for (layers in c("EFGHJP","E_G__P","_FGHJP","EFGH__")){
  for (envmode in c(0,1)){ 
    modelname <- switch(layers, "EFGHJP"="Full", "_FGHJP"="NoCue", "E_G__P"="NoHier", "EFGH__"="NoDev", "__G___"="Null"  )
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
boxatvec <- seq(1,12)
boxatvec <- boxatvec[-3*seq(1,4)]
axisatvec <- 3*seq(1,4)-1.5

#tiff("derr_AR.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
tiff(paste("ali_",pert,".tif",sep=""),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.ali100[2:ncol(df.ali100)], ylab = sprintf("Alignment (%s)",pert_vs), col=cols,ylim=c(0,1), at=boxatvec, xaxt = "n", cex.lab=1.5, cex.main=2.0)
axis(side = 1, at = axisatvec, labels=c("Full","NoHier","NoCue","NoDev"))
legend("topright", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
#legend("bottomleft",legend = c("Full","NoHier","NoCue","NoDev"), col=c("orange","navy","limegreen","darkorchid"), lty=c(rep(1,4)))
dev.off()

tiff(paste("sv1_",pert,".tif",sep=""),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.sv1100[2:ncol(df.sv1100)], ylab = sprintf("SV1 (%s)",pert_vs), col=colvec, at= boxatvec, xaxt="n", cex.lab=1.5, cex.main=2.0)
axis(side = 1, at = axisatvec, labels=c("Full","NoHier","NoCue","NoDev"))
legend("topright", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev"), col=c("orange","navy","limegreen","darkorchid"), lty=c(rep(1,4)))
dev.off()

tiff(paste("psv1_",pert,".tif",sep=""),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.psv1100[2:ncol(df.psv1100)], ylab = sprintf("%% 1st singular value (%s)",pert_vs), col=colvec, ylim=c(0,1), at = boxatvec, xaxt = "n", yaxt = "n", cex.lab=1.5, cex.main=2.0)

#To do: Fake up percentage axis values.

axis(side = 1, at = axisatvec, labels=c("Full","NoHier","NoCue","NoDev"))
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1.0), labels= c(0,20,40,60,80,100))
legend("topright", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
#legend("topright",legend = c("Full","NoHier","NoCue","NoDev"), col=c("orange","navy","limegreen","darkorchid"), lty=c(rep(1,4)))
dev.off()

tiff(paste("Fnorm_",pert,".tif",sep=""),width=2250,height=2250,units="px",pointsize=12,res=300)
#boxplot(df.Fnorm[2:5], main=paste(envmodename, pertname, sep=","), xlab="Model", ylab="Square Frobenius Norm", col=c("orange","limegreen","navy","darkorchid"), ylim=c(0,2.5))
boxplot(df.Fnorm[2:ncol(df.Fnorm)], ylab=sprintf("Total Cross Covariance (%s)",pert_vs), col=cols, at = boxatvec, xaxt = "n", cex.lab=1.5, cex.main=2.0) #May need to adjust range manually
axis(side = 1, at = axisatvec, labels=c("Full","NoHier","NoCue","NoDev"))
legend("topright", title="Environment", legend=c("Ancestral", "Novel"), lty=1, col=cols)
dev.off()