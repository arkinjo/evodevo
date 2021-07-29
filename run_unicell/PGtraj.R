
trajs <- read.table("traj1.dat",header=TRUE)
fit <- trajs$Fitness
dfit = diff(fit)
plot(dfit,main="Change in fitness",xlab="Generation",type="l")
abline(h=0,col="red")

nancs <- read.table("gen1_nanc.dat",header=TRUE)
nanc <- nancs$Ancestors

plot(nanc,xlab="Generation",ylim=c(0,1000),main="Shape of graph",type="l")


phat <- c()
ghat <- c()

pdf("PGtraj1.pdf")
for (i in c(1:200)){
  filename <- paste("PGtraj1_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  p <- pg$Phenotype
  g <- pg$Genotype
  plot(g,p,xlim=c(-15,15),ylim=c(-2,2),xlab="genotype",ylab="phenotype")
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

