
aphat <- c()
nphat <- c()
ghat <- c()

pdf("m5pop20210930hienvpg.pdf")
for (i in c(1:200)){
  filename <- paste("m5pop20210930hienvpg_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  ap <- pg$AncPhen
  np <- pg$NovPhen
  g <- pg$Genotype
  title <- paste("Generation",i,sep=" ")
  plot(g,ap,xlim=c(-10,10),ylim=c(-2,2),xlab="genotype",ylab="phenotype",main=title,col="red")
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