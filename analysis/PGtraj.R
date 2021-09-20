
aphat <- c()
nphat <- c()
ghat <- c()

pdf("upoppgtraj20210917-5.pdf")
for (i in c(1:200)){
  filename <- paste("upop20210917-5pg_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  ap <- pg$AncPhen
  np <- pg$NovPhen
  g <- pg$Genotype
  title <- paste("Generation",i,sep=" ")
  plot(g,ap,xlim=c(-5,5),ylim=c(-5,5),xlab="genotype",ylab="phenotype",main=title,col="red")
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