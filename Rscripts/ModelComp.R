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

df.ali <- data.frame(id=c(1:20))
df.sv1 <- data.frame(id=c(1:20))
df.psv1 <- data.frame(id=c(1:20))
df.Fnorm <- data.frame(id=c(1:20))

#df.play <- data.frame(id=c(1:20)) #Effective variance in direction of environmental change.

envmode <- 1 #Only consider novel environments

#denvs <- seq(10,100,10)
denvs <- append(2,seq(10,100,10))
#Challenge: Try to plot all models in one plot.
for (denv in denvs){
  for (layers in c("EFGHJP","_FGHJP")){ #Full, NoHier, NoCue, NoDev
    modelname <- switch(layers, "EFGHJP"="Full", "_FGHJP"="NoCue", "E_G__P"="NoHier", "EFGH__"="NoDev", "__G___"="Null"  )
    #modelenv <- paste(modelname,envmode)
    #envmodename <- switch(envmode+1, "Ancestral", "Novel", "Novel-Ancestral") 
    
    #Read data tables
    ali_in <- read.table(sprintf("%s_1_%s_ali.dat",layers,pert),header=FALSE)
    sv1_in <- read.table(sprintf("%s_1_%s_sval1.dat",layers,pert),header=FALSE)
    Fnorm_in <- read.table(sprintf("%s_1_%s_denv22.dat",layers,pert),header=FALSE)
    
    alivec <- ali_in$V2[which(sv1_in$V1==denv)]
    sv1vec <- sv1_in$V2[which(sv1_in$V1==denv)]

    modeldenv <- sprintf("%s_%d",modelname,denv)
    
    df.ali[modeldenv] <- alivec
    df.sv1[modeldenv] <- sv1vec
    df.psv1[modeldenv] <- sv1_in$V3[which(sv1_in$V1==denv)]
    df.Fnorm[modeldenv] <- Fnorm_in$V2[which(sv1_in$V1==denv)]
    
    #df.play[modeldenv] <- alivec*sv1vec

  }
  #df.ali100[modelenv] <- ali_in$V2[which(ali_in$V1==100)] #Alignment between first singular vector with environmental change
    
  #df.sv1100[modelenv] <- sv1_in$V2[which(sv1_in$V1==100)] #First singular value
  #df.psv1100[modelenv] <- sv1_in$V3[which(sv1_in$V1==100)] #Proportion first singular value
    
  #Fnorm_in <- read.table(paste(layers,"_",envmode,"_",pert,"_denv22.dat",sep=""),header=FALSE)
  #df.Fnorm[modelenv] <- Fnorm_in$V2[which(sv1_in$V1==100)] #Frobenius norm.
     
}


## Colors now indicate ancestral or novel environment
#cols <- c("magenta","cyan")
#colvec <- rep.int(cols,times=4)
#boxatvec <- seq(1,12)
#boxatvec <- boxatvec[-3*seq(1,4)]
#axisatvec <- 3*seq(1,4)-1.5

colvec <- c("red","cyan")
boxatvec <- seq(1,3*length(denvs))
boxatvec <- boxatvec[-3*seq(1,length(denvs))]
axisatvec <- seq(1.5,31.5,3.0)
pdenvs <- as.integer(denvs/2) #% change of environmental factors
spdenvs <- sprintf("%d",pdenvs) #string as boxplot label #Maybe not needed


png(sprintf("ali_%s.png",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.ali[2:ncol(df.ali)],ylab = sprintf("Alignment (%s)",pert_vs), xlab="% Environmental change",
        col=colvec, ylim=c(0,1), at=boxatvec, xaxt="n", cex.lab=1.5, cex.main=2.0)
axis(side=1, at=axisatvec, labels=spdenvs)
legend("topleft", title="Model", legend=c("Full","NoCue"), border="black", fill=colvec)
dev.off()

png(sprintf("sv1_%s.png",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
#make main plot first
par(fig=c(0,1,0,1)) #Specify NDC coordinates for main plot
boxplot(df.sv1[2:ncol(df.sv1)],ylab = sprintf("1st singular value (%s)",pert_vs), xlab="% Environmental change", 
        col=colvec, at=boxatvec, xaxt="n", cex.lab=1.5, cex.main=2.0)
axis(side=1, at=axisatvec, labels=spdenvs)
legend("topleft", title="Model", legend=c("Full","NoCue"), border="black", fill=colvec)
#make inset plot
par(fig=c(0.05,0.65,0.3,0.8),new=T)
boxplot(df.sv1[2:13], col=colvec, at=boxatvec[1:12], xaxt="n", log="y", cex.axis=0.8) #log scale
axis(side=1, at=axisatvec, labels=spdenvs, cex.axis=0.8)
dev.off()

png(sprintf("psv1_%s.png",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
boxplot(df.psv1[2:ncol(df.psv1)],ylab = sprintf("%% 1st singular value (%s)",pert_vs), xlab="% Environmental change",
        col=colvec, ylim=c(0,1), at=boxatvec, xaxt="n", yaxt="n", cex.lab=1.5, cex.main=2.0)
axis(side = 1, at = axisatvec, labels = spdenvs)
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1.0), labels= c(0,20,40,60,80,100))
legend("topleft", title="Model", legend=c("Full","NoCue"), border="black", fill=colvec)
dev.off()

png(sprintf("Fnorm_%s.png",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
par(fig=c(0,1,0,1)) #Specify NDC coordinates for main plot
boxplot(df.Fnorm[2:ncol(df.Fnorm)],ylab = sprintf("1st singular value (%s)",pert_vs), xlab="% Environmental change",
        col=colvec, at=boxatvec, xaxt="n", cex.lab=1.5, cex.main=2.0)
axis(side=1, at=axisatvec, labels=spdenvs)
legend("topleft", title="Model", legend=c("Full","NoCue"), border="black", fill=colvec)
#make inset plot
par(fig=c(0.05,0.65,0.3,0.8),new=T)
boxplot(df.Fnorm[2:13], col=colvec, at=boxatvec[1:12], xaxt="n", log="y", cex.axis=0.8) #log scale
axis(side=1, at=axisatvec, labels=spdenvs, cex.axis=0.8)

dev.off()

#png(sprintf("play_%s.png",pert),width=2250,height=2250,units="px",pointsize=12,res=300)
#boxplot(df.play[2:ncol(df.play)],ylab = sprintf("1st singular value (%s)",pert_vs), xlab="% Environmental change",
#        col=colvec, ylim=c(0,1), at=boxatvec, xaxt="n", cex.lab=1.5, cex.main=2.0)
#axis(side=1, at=axisatvec, labels=spdenvs)
#legend("topright", title="Model", legend=c("Full","NoHier","NoCue","NoDev"), lty=1, col=cols)
#dev.off()

