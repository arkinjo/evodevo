
#Run this after running updated gvar.R and traj.R scripts!

df.err_minvar <- data.frame(epoch=seq(1,nepoch)) #error at minimum genetic variance
df.derr_GA <- data.frame(epoch=seq(1,nepoch)) #change in error due to change in regulation
df.derr_AR <- data.frame(epoch=seq(1,nepoch)) #change in error due to adaptive refinement

df.pderr_GA <- data.frame(epoch=seq(1,nepoch))
df.pderr_AR <- data.frame(epoch=seq(1,nepoch)) #Proportion change in error due to adaptive refinement. 1- gives proportion due to change in regulation

df.derr_ALL <- data.frame(epoch=seq(1,nepoch))

for(i in seq(2,ncol(df.tminvar))){
  err_minvar <- c()
  derr_GA <- c()
  derr_AR <- c()
  pderr_AR <- c()
  model <- colnames(df.tminvar)[i]
  modelgens <- df.tminvar[model]
  for(j in seq(1,nrow(df.tminvar))){
    modelXepoch <- sprintf("%s_%02d",model,j)
    gen_minvar <- modelgens[j,]
    errtraj <- df.err[modelXepoch]
    err_minvar <- append(err_minvar,errtraj[gen_minvar,])
    derr_GA <- append(derr_GA,errtraj[1,]-errtraj[gen_minvar,])
    derr_AR <- append(derr_AR,errtraj[gen_minvar,]-errtraj[ngen,])
    pderr_AR <- append(pderr_AR,(errtraj[gen_minvar,]-errtraj[ngen,])/(errtraj[1,]-errtraj[ngen,]))
  }
  df.err_minvar[model] <- err_minvar
  df.derr_GA[model] <- derr_GA
  df.derr_ALL[sprintf("%s_GA",model)] <- derr_GA
  df.derr_AR[model] <- derr_AR
  df.derr_ALL[sprintf("%s_AR",model)] <- derr_AR
  df.pderr_AR[model] <- pderr_AR
  df.pderr_GA[model] <- 1-pderr_AR
}

#Make boxplots of change in error in each phase

boxatvec <- seq(1,12)
boxatvec <- boxatvec[-3*seq(1,4)]
axisatvec <- 3*seq(1,4)-1.5

#No need length of adaptive refinement phase, as long as environment doesn't change this could be infinite

tiff("derr_GA.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.derr_GA[bxpindex+1],col=colvec.s[bxpindex],ylab="Decrease in mismatch (Change in regulation)",cex.lab=1.5,cex.main=2.0)
dev.off()

tiff("derr_AR.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.derr_AR[bxpindex+1],col=colvec.s[bxpindex],ylab="Decrease in mismatch (Adaptive refinement)",cex.lab=1.5,cex.main=2.0)
dev.off()

tiff("pderr_AR.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.pderr_AR[bxpindex+1],col=colvec.s[bxpindex],ylab="Proportion decrease in mismatch",main="Adaptive refinement")
dev.off()

tiff("pderr_GA.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.pderr_GA[bxpindex+1],col=colvec.s[bxpindex],xlab="Model",ylab="Proportion decrease in error",main="Change in regulation")
dev.off()

tiff("derr_ALL.tif",width=2250,height=2250,units="px", pointsize=12, res=300)
boxplot(df.derr_ALL[c(6,7,4,5,2,3,8,9)],ylim=c(0,1),col=c("magenta","cyan"),ylab="Decrease in mismatch",at=boxatvec,xaxt="n",cex.lab=1.5)
axis(side = 1, at = axisatvec, labels=c("Full","NoHier","NoCue","NoDev"))
legend("topright",legend=c("Before bottleneck","After bottleneck"),col=c("magenta","cyan"),lty=1)
dev.off()

# #For fun
# df.tminvar1 <- df.tminvar[-1]
# tminvarvec <- unlist(df.tminvar1)
# durvec <- ngen - tminvarvec # Number of generations in adaptive refinement phase
# df.err_inf1 <- df.err_inf[-1]
# err_infvec <- unlist(df.err_inf1)
# 
# df.tVSerr <- data.frame(t=tminvarvec, dur=durvec, err=err_infvec)
# lm.tVSerr <- lm(df.tVSerr$err~df.tVSerr$t) #Linear regression
# summary(lm.tVSerr) #Slope insignificant.
# 
# 
# png("play.png")
# plot(x=df.tVSerr$t,y=df.tVSerr$err,col=colvec.l,pch=c(rep(1,nepoch),rep(2,nepoch),rep(3,nepoch),rep(4,nepoch)),xlab="Generations to bottleneck",ylab="Final Error")
# legend("topright",legend=c("Full","NoHier","NoCue","NoDev"), pch=bxpindex, col=colvec.s[bxpindex])
# dev.off()
# 
# png("play2.png")
# plot(x=df.tVSerr$dur,y=df.tVSerr$err,col=colvec.l,pch=c(rep(1,nepoch),rep(2,nepoch),rep(3,nepoch),rep(4,nepoch)),xlab="Adaptive refinement duration",ylab="Final Error")
# legend("topleft",legend=c("Full","NoHier","NoCue","NoDev"), pch=bxpindex, col=colvec.s[bxpindex])
# dev.off()