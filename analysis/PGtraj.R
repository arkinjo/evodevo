
phat <- c()
ghat <- c()

pdf("upoppgtraj20210913.pdf")
for (i in c(0:201)){
  filename <- paste("upop20210913pg_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  p <- pg$Phenotype
  g <- pg$Genotype
  title <- paste("Generation",i,sep=" ")
  plot(g,p,xlim=c(-15,15),ylim=c(-50,50),xlab="genotype",ylab="phenotype",main=title)
  abline(h=0)
  abline(v=0)
  mup = mean(p)
  mug = mean(g)
  #print(mup)
  #print(mug)
  phat = append(phat,mup)
  ghat = append(ghat,mug)
}
dev.off()

dp = diff(phat)
dg = diff(ghat)

plot(dp,main="Change in phenotype", xlab="Generation",type="l")
abline(h=0,col="red")
plot(dg,main="Change in genotype", xlab="Generation",type="l")
abline(h=0,col="red")

