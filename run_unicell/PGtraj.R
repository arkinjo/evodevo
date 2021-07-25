
trajs <- read.table("traj4.dat",header=TRUE)
fit <- trajs$Fitness
dfit = diff(fit)
plot(dfit,main="Change in fitness",xlab="Generation",type="l")
abline(h=0,col="red")

phat <- c()
ghat <- c()

for (i in c(1:200)){
  filename <- paste("PGtraj4_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  p <- pg$Phenotype
  g <- pg$Genotype
  pdfout <- paste("PGtraj4_",i,".pdf",sep="")
  pdf(pdfout)
  plot(g,p,xlim=c(-15,15),ylim=c(-2,2),xlab="genotype",ylab="phenotype")
  dev.off()
  phat = append(phat,mean(p))
  ghat = append(ghat,mean(g))
}

dp = abs(diff(phat))
dg = abs(diff(ghat))

plot(dp,main="Change in phenotype", xlab="Generation",type="l")
plot(dg,main="Change in genotype", xlab="Generation",type="l")