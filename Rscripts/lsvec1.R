
#envmode <- 0 ### 0: Ancestral, 1: Novel, 2: Diff
#pert <- "pe" ### "pG": Genetic/Mutational perturbation, "pe": environmental perturbation

#envmodename <- switch(envmode+1, "Ancestral", "Novel", "Novel-Ancestral") 
#pertname <- switch(pert, "pG" = "Mutation", "pe" = "Environmental Change") 

df.eGali <- data.frame(id=c(1:20))
colvec <- c()

for (layers in c("EFGHJP","E_G__P","_FGHJP","EFGH__")){
  modelname <- switch(layers, "EFGHJP"="Full", "_FGHJP"="NoCue", "E_G__P"="NoHier", "EFGH__"="NoDev", "__G___"="Null"  )
  for(envmode in c(0,1)){
    modelXenvmode <- sprintf("%s_%s",modelname,envmode)
    envcolor <- switch(envmode+1,"magenta","cyan")
    eGalilist <- c() #Initialize for each model
    for (epoch in c(1:20)){
      str.epoch <- sprintf("%02d",epoch)
      filename <- paste(layers,"_100_0.05_",envmode,"_",str.epoch,".lsvec1",sep="")
      df.lsvec <- read.table(filename,header=FALSE)
      evec <- df.lsvec$V2
      Gvec <- df.lsvec$V3
      eGalilist <- append(abs(sum(evec*Gvec)),eGalilist) #Calculate absolute value of dot product
    }
    df.eGali[modelXenvmode] <- eGalilist
    colvec <- append(colvec,envcolor)
  }
 
}

boxatvec <- c(1:12) 
boxatvec <- boxatvec[-3*c(1:4)]
axisatvec <- c(1.5,4.5,7.5,10.5)
png("eGali.png")
boxplot(df.eGali[2:ncol(df.eGali)], ylab="Pheno(Cue)-Pheno(Geno) alignment (Gen 1)", ylim=c(0,1),col=colvec,xaxt="n",at=boxatvec,cex.lab=1.5)
axis(1,at=axisatvec,labels=c("Full","NoHier","NoCue","NoDev"))
legend("bottomleft",legend=c("Ancestral","Novel"),col = c("magenta","cyan"),lty=1, title = "Environment")
dev.off()

#Plot both in same environment
