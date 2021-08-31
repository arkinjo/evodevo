
phat <- c()
ghat <- c()

pdf("m2poptraj20210829pg.pdf")
for (i in c(0:200)){
  filename <- paste("m2pop20210829pg_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  p <- pg$Phenotype
  g <- pg$Genotype
  title <- paste("Generation",i,sep=" ")
  plot(g,p,xlim=c(-15,15),ylim=c(-3,3),xlab="genotype",ylab="phenotype",main=title)
  abline(h=0)
  abline(v=0)
  phat = append(phat,mean(p))
  ghat = append(ghat,mean(g))
}
dev.off()

dp = diff(phat)
dg = diff(ghat)

plot(dp,main="Change in phenotype", xlab="Generation",type="l")
abline(h=0,col="red")
plot(dg,main="Change in genotype", xlab="Generation",type="l")
abline(h=0,col="red")

