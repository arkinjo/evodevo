
ptraj <- read.table("phen1_merged.dat", header=FALSE)
c.ptraj <- scale(ptraj,scale=FALSE)
mu <- attr(c.ptraj,"scaled:center")

svd.ptraj <- svd(c.ptraj)

xaxis <- svd.ptraj$v[,1]
yaxis <- svd.ptraj$v[,2]

for (gen in c(1:200)){
  filename = paste("phen1_",gen,".dat",sep="")
  phen <- read.table(filename,header=FALSE)
  c.phen = phen-mu 
  mat.c.phen = as.matrix(c.phen)
  x = c(1:length(mat.c.phen[,1]))
  for (i in c(1:length(mat.c.phen[,1]))){
    x[i] = sum(xaxis*mat.c.phen[i,])
  }
  y = c(1:length(mat.c.phen[,1]))
  for (i in c(1:length(mat.c.phen[,1]))){
    y[i] = sum(yaxis*mat.c.phen[i,])
  }
  title = paste("Generation",gen,sep="_")
  plot(x,y,type="p",xlab="PC1",ylab="PC2",xlim=c(-3,3),ylim=c(-3,3),main=title)
  abline(h=0)
  abline(v=0)
}

svs <- svd.ptraj$d
plot(svs,main="Plot of square of singular values",xlab="Index",ylab="Singular value")
plot(cumsum(svs),main="Cumulative sum of squares",xlab="Index",ylab="Sum")

svs2 <- svs*svs

plot(svs2,main="Plot of square of singular values",xlab="Index",ylab="Singular value")
plot(cumsum(svs2),main="Cumulative sum of squares",xlab="Index",ylab="Sum")
