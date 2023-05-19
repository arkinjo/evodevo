### Training data

traj <- read.table("EFGHJP_train.traj",header=FALSE)

nepoch <- max(traj$V1)


for (epoch in c(1:nepoch)) {
  row_index <- which(traj$V1 == epoch) #Extract rows
  
  gen <- traj$V2[row_index]
  errtraj <- traj$V5[row_index]
  
  tiff(sprintf("err_train%02d.tif",epoch),width=2250,height=2250,units="px", pointsize=12, res=300)
  plot(gen,errtraj,xlab="Generation",ylab="Mismatch",main=sprintf("Training Epoch %d",epoch),type="l",cex.lab=1.5,cex.main=2.0)
  dev.off()
}