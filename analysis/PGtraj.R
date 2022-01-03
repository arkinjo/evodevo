
aphat <- c()
nphat <- c()
ghat <- c()

pdf("upop20211225-30testpg.pdf")
for (i in c(1:200)){
  filename <- paste("upop20211225-30testpg_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  ap <- pg$AncPhen
  np <- pg$NovPhen
  g <- pg$Genotype
  title <- paste("Generation",i,sep=" ")
  plot(g,ap,xlim=c(-10,10),ylim=c(-10,10),xlab="genotype",ylab="phenotype",main=title,col="red")
  points(g,np,col="blue")
  abline(h=0)
  abline(v=0)
  muap = mean(ap)
  munp = mean(np)
  mug = mean(g)
  aphat = append(aphat,muap)
  nphat = append(nphat,munp)
  ghat = append(ghat,mug)
}
dev.off()

plot(nphat,ylim=c(-10,10),type="l",xlab="Generation",ylab="Projected Phenotype",main="Mean phenotype",col="blue")
abline(h=0)
points(aphat,type="l",col="red")

