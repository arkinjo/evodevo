### Plot phenotype expression against genotype
g0 <- read.table("gen1_1.dat",header=FALSE)
cg0 <- scale(g0,scale=FALSE)
mu_g0 <- attr(cg0,"scaled:center")

g1 <- read.table("gen1_200.dat",header=FALSE)
cg1 <- scale(g1,scale=FALSE)
mu_g1 <- attr(cg1,"scaled:center")

haxis = mu_g1-mu_g0
eh = haxis/sqrt(haxis*haxis)

ptraj <- read.table("phen1_merged.dat",header=FALSE)
cptraj <- scale(ptraj,scale=FALSE)
mu <- attr(cptraj,"scaled:center")
svd.ptraj <- svd(cptraj)
vaxis = svd.ptraj$v[,1]

#To do: get centre of phenotype and genotype trajectory

xbar = c()
ybar = c()

for (gen in c(1:200)){
  gfilename <- paste("gen1_",gen,".dat",sep="")
  genomes <- read.table(gfilename,header=FALSE) #centralizing not needed for genome
  mat.cgen <- as.matrix(genomes)
  x <- c(1:length(mat.cgen[,1]))
  for (i in c(1:length(mat.cgen[,1]))){
    x[i] = sum(haxis*mat.cgen[i,])
  }
  
  pfilename <- paste("phen1_",gen,".dat",sep="")
  phen <- read.table(pfilename,header=FALSE)
  cphen = phen-mu
  mat.cphen <- as.matrix(cphen)
  y <- c(1:length(mat.cphen[,1]))
  for (i in c(1:length(mat.cphen[,1]))){
    y[i] = sum(vaxis*mat.cphen[i,])
  }
  
  title = paste("Generation",gen,sep=" ")
  pdftitle = paste("PGtraj1_",gen,".pdf",sep="")
  pdf(pdftitle)
  plot(x,y,main=title,xlab="genotype",ylab="phenotype",xlim=c(-300,300),ylim=c(-3,3))
  abline(h=0)
  abline(v=0)
  dev.off()
  xbar = append(xbar,mean(x))
  ybar = append(ybar,mean(y))
}

pdf("meanPGtraj")
plot(xbar,ybar,main="mean trajectory",xlab="genotype",ylab="phenotype",xlim=c(-300,300),ylim=c(-3,3))
dev.off()

dx = abs(diff(xbar)) #Magnitude of change
dy = abs(diff(ybar)) 
pdf("diffGtraj")
plot(dx,main="Change in genotype",xlab="Generation",ylab="Magnitdue of change")
dev.off()
pdf("diffPtraj")
plot(dy,main="Change in phenotype",xlab="Generation",ylab="Magnitdue of change")
dev.off()
