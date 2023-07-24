
#layerlist <- c("EFGHJP","E_G__P","_FGHJP")

#Make sure to set working directory to source file location prior to 
setwd("../proj")

#par(cex=1.0, cex.axis=1.0, cex.lab=1.7, cex.main=2.0) #Warning! This also carries forward to legend font size


maxgen <- 200 #Epoch length
nepoch <- 10
#fpt.Threshold <- 0.9 #Arbitrary number between 0 and 1.
#str.fpt <- sprintf("%.2f",fpt.Threshold)

df.G <- data.frame(gen=c(1:maxgen))
df.dG <- data.frame(gen=c(1:(maxgen-1)))
npvar <- data.frame(gen=c(1:maxgen))
apvar <- data.frame(gen=c(1:maxgen))

df.pGvar <- data.frame(gen=c(1:maxgen))

df.p_0 <- data.frame(epoch=seq(1,nepoch)) #Projected phenotypes before evolution
df.p_inf <- data.frame(epoch=seq(1,nepoch)) #Projected phenotypes after evolution


#models <- c("Full","NoCue","NoHier","NoDev")
modelvec <- c()
colvec.s <- c()
colvec.l <- c()

for(layers in c("_FGHJP","E_G__P","EFGHJP","EFGH__","CEFGHJP")){ #Loop over models, order determines what goes over/under in plot! #Ignore NoTrain in main plots.
  
  modelname <- switch(layers, "EFGHJP"="Full","_FGHJP"="NoCue","E_G__P"="NoHier","EFGH__"="NoDev","CEFGHJP"="NoTrain")
  modelvec <- append(modelvec, modelname)
  modelcol <- switch(layers, "EFGHJP" = "orange", "_FGHJP"="cyan", "E_G__P"="limegreen", "EFGH__"="darkorchid", "CEFGHJP"="blue")
  colvec.s <- append(colvec.s, modelcol)
  
  ap_0 <- c()
  np_0 <- c()
  ap_inf <- c()
  np_inf <- c()
  
  for(pertsize in c(10)) { #Genotype-phenotype plot in small and large environmental changes
    for(epoch in c(1:nepoch)) {
      colvec.l <- append(colvec.l,modelcol)
    
      aghat <- c()
      aphat <- c()
      nghat <- c()
      nphat <- c()
    
      gvar <- c()
    
      sdap <- c()
      sdnp <- c()
      sdag <- c()
      sdng <- c()
    
      print(sprintf("Plotting ../anims/%s_%d_%02d.pdf",modelname,pertsize,epoch))
      pdf(sprintf("../anims/%s_%d_%02d.pdf",modelname,pertsize,epoch)) #dump animation in seperate directory
      for (gen in c(1:maxgen)){
        t0gen <- sprintf("%03d",gen)
        #filename <- paste(layers,"_run100_",str.epoch,"_",t0gen,".dat",sep="")
        filename <- sprintf("%s_run%d_%02d_%03d.dat",layers,pertsize,epoch,gen)
      
        pgstats <- read.delim(filename,skip = 1000)
        #print(filename)
      
        aghat <- append(aghat,pgstats[1,3])
        aphat <- append(aphat,pgstats[1,4])
        sdag <- append(sdag,pgstats[1,5])
        sdap <- append(sdap,pgstats[1,6])
        nghat <- append(nghat,pgstats[2,3])
        nphat <- append(nphat,pgstats[2,4])
        sdng <- append(sdng,pgstats[2,5])
        sdnp <- append(sdnp,pgstats[2,6])
      
        pgpoints <- read.delim(filename,nrows = 1000)
        plot(pgpoints$Geno.e0,pgpoints$Pheno0,xlim=c(0.0,1.0),ylim=c(0.0,1.0),
             xlab="Genotype",ylab="Phenotype",main=sprintf("%6s (Gen:%3d)", modelname, gen),
             col="darkorchid",pch=1,cex.lab=1.5,cex.main=2.0)
        points(pgpoints$Geno.e1,pgpoints$Pheno1,col="cyan",pch=4)
        #legend("topleft",legend=c("Novel","Ancestral"),col=c("cyan","darkorchid"),pch=c(4,1), title="Environment") 
        #Can't find a good place to put figure legend without getting into way of plot
      }
      dev.off()
    
      #Summarize trajectory: error plots for initial, reference and final generation.
      #pdf(paste("../plots/",modelname,"_",str.epoch,".pdf",sep=""))
    
    
      png(sprintf("../plots/%s_%d_%02d.png",modelname,pertsize,epoch),width=2250,height=2250,units="px", pointsize=12, res=300)
      plot(aghat,aphat,xlim=c(0.0,1.0),ylim=c(0.0,1.0),xlab="Genotype",ylab="Phenotype",main=modelname,col="darkorchid",type="l",cex.lab=1.5,cex.main=2.0)
      points(aghat,aphat,col="darkorchid",pch=1)
      arrows(aghat, aphat-sdap, aghat, aphat+sdap, length=0.05, angle=90, code=3, col="darkorchid") 
      arrows(aghat-sdag, aphat, aghat+sdag, aphat, length=0.05, angle=90, code=3, col="darkorchid") 
      lines(nghat,nphat,col="cyan")
      points(nghat,nphat,col="cyan",pch=4)
      arrows(nghat, nphat-sdnp, nghat, nphat+sdnp, length=0.05, angle=90, code=3, col="cyan") 
      arrows(nghat-sdng, nphat, nghat+sdng, nphat, length=0.05, angle=90, code=3, col="cyan") 
      legend("topleft",legend=c("Novel","Ancestral"),col=c("cyan","darkorchid"),pch=c(4,1), title="Environment",lty=c(1,1))
      dev.off()
    
      #df.G[paste(modelname,str.epoch,sep="")] <- nghat
      #df.dG[paste(modelname,str.epoch,sep="")] <- diff(nghat)
      #apvar[paste(modelname,str.epoch,sep="")] <- sdap^2
      #npvar[paste(modelname,str.epoch,sep="")] <- sdnp^2
    
      #df.pGvar[sprintf("%s%02d",modelname,epoch)] <- sdng^2
    
      #ap_0 <- append(ap_0,aphat[1])
      #np_0 <- append(np_0,nphat[1])
      #ap_inf <- append(ap_inf,aphat[maxgen])
      #np_inf <- append(np_inf,nphat[maxgen])
    
    
      print(sprintf("Plotting ../plots/%s_%d_%02d.png",modelname,pertsize,epoch))
    }

  }
  
  #df.p_0[sprintf("%s_a",modelname)] <- ap_0
  #df.p_0[sprintf("%s_n",modelname)] <- np_0
  #df.p_inf[sprintf("%s_a",modelname)] <- ap_inf
  #df.p_inf[sprintf("%s_n",modelname)] <- np_inf
  
  #df.fpt.ag[modelname] <- fpt.ag
  #df.fpt.ap[modelname] <- fpt.ap
  #df.fpt.ng[modelname] <- fpt.ng
  #df.fpt.np[modelname] <- fpt.np
  
}

#Change working directory to NoDev here! #No need anymore
#setwd("../nodev")

mat.G <- as.matrix(df.G)
#dG <- diff(df.G) #Bad line
mat.dG <- as.matrix(df.dG)

mat.apvar <- as.matrix(apvar)
mat.npvar <- as.matrix(npvar)

mat.pGvar <- as.matrix(df.pGvar)

lindex <- c(3,2,1,4)
boxatvec <- seq(1,12)
boxatvec <- boxatvec[-3*seq(1,4)]
axisatvec <- 3*seq(1,4)-1.5


