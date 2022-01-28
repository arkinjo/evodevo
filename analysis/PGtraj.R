
refgen <- 20 #Reference generation
maxgen <- 200 #Epoch length
popseed <- 20220115
refgens <- c(1,refgen,maxgen)

aphat <- c()
nphat <- c()
ghat <- c()
sdap <- c()
sdnp <- c()
sdg <- c()

pdffilename <- paste("upop",popseed,"testpg.pdf",sep="")
pdf(pdffilename)
for (i in c(1:maxgen)){
  filename <- paste("upop",popseed,"testpg_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  ap <- pg$AncPhen
  np <- pg$NovPhen
  g <- pg$Genotype
  title <- paste("Generation",i,sep=" ")
  plot(g,ap,xlim=c(-20,20),ylim=c(-10,10),xlab="genotype",ylab="phenotype",main=title,col="red")
  points(g,np,col="blue")
  abline(h=0)
  abline(v=0)
  #muap = mean(ap)
  #munp = mean(np)
  #mug = mean(g)
  aphat = append(aphat,mean(ap))
  nphat = append(nphat,mean(np))
  ghat = append(ghat,mean(g))
  sdap = append(sdap,sd(ap))
  sdnp = append(sdnp,sd(np))
  sdg = append(sdg,sd(g))
}
dev.off()

#Summarize trajectory: error plots for initial, reference and final generation.
plot(ghat,aphat,xlim=c(-20,20),ylim=c(-10,10),xlab="Genotype",ylab="Phenotype",main="Mean Projected Trajectory",col="red",type="l")
points(ghat,aphat,col="red")
arrows(ghat[refgens], aphat[refgens]-sdap[refgens], ghat[refgens], aphat[refgens]+sdap[refgens], length=0.05, angle=90, code=3, col="red") 
arrows(ghat[refgens]-sdg[refgens], aphat[refgens], ghat[refgens]+sdg[refgens], aphat[refgens], length=0.05, angle=90, code=3, col="red") 
lines(ghat,nphat,col="blue")
points(ghat,nphat,col="blue")
arrows(ghat[refgens], nphat[refgens]-sdnp[refgens], ghat[refgens], nphat[refgens]+sdnp[refgens], length=0.05, angle=90, code=3, col="blue") 
arrows(ghat[refgens]-sdg[refgens], nphat[refgens], ghat[refgens]+sdg[refgens], nphat[refgens], length=0.05, angle=90, code=3, col="blue") 


