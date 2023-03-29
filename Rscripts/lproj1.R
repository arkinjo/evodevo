
#envmode <- 0 ### 0: Ancestral, 1: Novel, 2: Diff
#pert <- "pe" ### "pG": Genetic/Mutational perturbation, "pe": environmental perturbation

#envmodename <- switch(envmode+1, "Ancestral", "Novel", "Novel-Ancestral") 
#pertname <- switch(pert, "pG" = "Mutation", "pe" = "Environmental Change") 

#df.ccaproj <- data.frame(id=c(1:1000))
colvec <- c()

for (layers in c("EFGHJP","E_G__P","_FGHJP")){
  modelname <- switch(layers, "EFGHJP"="Full", "_FGHJP"="NoCue", "E_G__P"="NoHier", "EFGH__"="NoDev", "__G___"="Null"  )
  
  for(epoch in c(1:20)){
    #modelXenvmode <- sprintf("%s_%s",modelname,envmode)
    #envcolor <- switch(envmode+1,"magenta","cyan")
    eGalilist <- c() #Initialize for each model
    #for (envmode in c(0,1)){
    #str.epoch <- sprintf("%02d",epoch)
    afilename <- sprintf("%s_100_0.05_0_%02d.lproj1",layers,epoch)
    df.alproj1 <- read.table(afilename,header=FALSE)
    eproj_a <- df.alproj1$V2
    Gproj_a <- df.alproj1$V3
    cor_a <- cor(eproj_a,Gproj_a)
    #if(cor_a < 0.0){ #Flip
    #  Gproj_a <- -Gproj_a
    #}
    
    afilename <- sprintf("%s_100_0.05_1_%02d.lproj1",layers,epoch)
    df.nlproj1 <- read.table(afilename,header=FALSE)
    eproj_n <- df.nlproj1$V2
    Gproj_n <- df.nlproj1$V3
    cor_n <- cor(eproj_n,Gproj_n)
    #if(cor_n < 0.0){ #Flip
    #  Gproj_n <- -Gproj_n
    #}
    
    tiff(sprintf("%s-ccaproj-%02d.tif",modelname,epoch),width=2250,height=2250,units="px",pointsize=12,res=300)
    plot(Gproj_n,eproj_n,col="cyan",xlab="Pheno-Geno projection",ylab="Pheno-Cue projection",main=modelname,
         cex.lab=1.5,cex.main=2.0,xlim=c(-10,10),ylim=c(-10,10))
    abline(h=0)
    abline(v=0)
    points(Gproj_a,eproj_a,col="magenta")
    legend("bottomright",legend=c(sprintf("Ancestral (cor = %.5f)",cor_a),sprintf("Novel (cor = %.5f)",cor_n)),
                                  lty=1,col=c("magenta","cyan"),title="Environment")
    dev.off()
    #eGalilist <- append(abs(sum(evec*Gvec)),eGalilist) #Calculate absolute value of dot product
    #df.eGali[modelXenvmode] <- eGalilist
    #colvec <- append(colvec,envcolor)
  }
  
}

#Need: autocorrect to positive correlation by multiplying by -1.

