#Do one model at a time first

pert <- "pG" ### "pG": Genetic/Mutational perturbation, "pe": environmental perturbation
layers <- "E_G__P" #Layers present in model
modelname <- switch(layers, "EFGHJP"="Full", "_FGHJP"="NoCue", "E_G__P"="NoHier", "EFGH__"="NoDev", "__G___"="Null")
#Reminder that R starts counting at 1 instead of 0.
pert_vs <- switch(pert, "pG" = "Pheno-Geno", "pe" = "Pheno-Cue") 

df.ali <- data.frame(id=c(1:20))
df.sv1 <- data.frame(id=c(1:20))
df.psv1 <- data.frame(id=c(1:20))
df.Fnorm <- data.frame(id=c(1:20))

envmode <- 1 #Only consider novel environments

#denvs <- seq(10,100,10)

#Read data tables

ali_in <- read.table(sprintf("%s_1_%s_ali.dat",layers,pert),header=FALSE)
sv1_in <- read.table(sprintf("%s_1_%s_sval1.dat",layers,pert),header=FALSE)
Fnorm_in <- read.table(sprintf("%s_1_%s_denv22.dat",layers,pert),header=FALSE)

denvs <- append(2,seq(10,100,10))
#Challenge: Try to plot all models in one plot. #Not here, do model by model first.
for (denv in denvs){
  
  #for (layers in c("EFGHJP","E_G__P","_FGHJP","EFGH__")){ #Full, NoHier, NoCue, NoDev
    
  pdenv <- as.integer(denv/2) #Percentage change in environmental factors
  spdenv <- sprintf("%d",pdenv) #string as boxplot label
  
  df.ali[spdenv] <- ali_in$V2[which(ali_in$V1==denv)]
  df.sv1[spdenv] <- sv1_in$V2[which(sv1_in$V1==denv)]
  df.psv1[spdenv] <- sv1_in$V3[which(sv1_in$V1==denv)]
  df.Fnorm[spdenv] <- Fnorm_in$V2[which(Fnorm_in$V1==denv)]
  
}

png(sprintf("ali_%s_%s.png",modelname,pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.ali[2:ncol(df.ali)],ylab = sprintf("Alignment (%s)",pert_vs), ylim=c(0,1), xlab="% Environmental change", 
        main=modelname, col="magenta", cex.lab=1.5, cex.main=2.0)
#axis()
#legend("topright", title="Model", legend=c("Full","NoHier","NoCue","NoDev"), lty=1, col=cols)
dev.off()

png(sprintf("sv1_%s_%s.png",modelname,pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.sv1[2:ncol(df.sv1)],ylab = sprintf("1st singular value (%s)",pert_vs), xlab="% Environmental change",
        main=modelname, col="magenta", cex.lab=1.5, cex.main=2.0)
#axis()
#legend("topright", title="Model", legend=c("Full","NoHier","NoCue","NoDev"), lty=1, col=cols)
dev.off()

png(sprintf("psv1_%s_%s.png",modelname,pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.psv1[2:ncol(df.psv1)],ylab = sprintf("%% 1st singular value (%s)",pert_vs), ylim=c(0,1), xlab="% Environmental change",
        main=modelname, yaxt="n", col="magenta", cex.lab=1.5, cex.main=2.0)
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1.0), labels= c(0,20,40,60,80,100))
#legend("topright", title="Model", legend=c("Full","NoHier","NoCue","NoDev"), lty=1, col=cols)
dev.off()

png(sprintf("Fnorm_%s_%s.png",modelname,pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.Fnorm[2:ncol(df.sv1)],ylab = sprintf("Total Cross-Covariance (%s)",pert_vs), xlab="% Environmental change",
        main=modelname, col="magenta", cex.lab=1.5, cex.main=2.0)
#axis()
#legend("topright", title="Model", legend=c("Full","NoHier","NoCue","NoDev"), lty=1, col=cols)
dev.off()

