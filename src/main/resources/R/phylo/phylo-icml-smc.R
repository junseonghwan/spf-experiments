rm(list=ls())
library(spatstat)

dir<-paste("~/Dropbox/Research/repo/spf-experiments/output/icml3/smc/", sep="")
folders1<-c("10", "20", "40", "50")
method_type<-"smc"
N<-length(folders1)
idx<-1
dd<-matrix(0, ncol=1, nrow = length(folders1))
hh<-matrix(0, ncol=1, nrow = length(folders1))
logZ<-matrix(0, ncol=1, nrow = length(folders1))
pdf("~/temp/icml-smc.pdf",width=6,height=4,paper='special')
for (i in 1:length(folders1))
{
  folder1<-folders1[i]
  d1k<-read.csv(paste(dir, folder1, "/output", idx, "/heights-", method_type, ".csv", sep=""), header=F)
  ww<-read.csv(paste(dir, folder1, "/output", idx, "/weights-", method_type, ".csv", sep=""), header=F)
  truth<-read.csv(paste(dir, folder1, "/output", idx, "/height-truth.csv", sep=""), header=F)
  hh[i,1]<-truth$V1
  
  plot(density(d1k$V1, weights = ww$V1), col='black', xlim=c(min(c(d1k$V1, truth$V1)), max(c(d1k$V1, truth$V1))), main = paste("Simulation", i))
  abline(v=truth$V1, col='red')
  temp<-ewcdf(d1k$V1, ww$V1)
  q1k<-quantile(temp, probs = c(0.025, 0.975))
  abline(v=q1k, col='black')

  trueDistMatrix<-read.csv(paste(dir, folder1, "/output", idx, "/pairwise-dist-smc.csv", sep=""), header=T)
  distMatrix<-read.csv(paste(dir, folder1, "/output", idx, "/pairwise-dist-truth.csv", sep=""), header=T)
  plot(hclust(dist(trueDistMatrix/2)), main="true tree")
  plot(hclust(dist(distMatrix/2)), main="mean tree")
  diff<-trueDistMatrix - distMatrix
  dd[i,1]<-sum(diff[upper.tri(diff)]^2)
  
  logZ[i,1]<-read.csv(paste(dir, folder1, "/output", idx, "/logZ-", method_type, ".csv", sep=""), header=F)$V1[1]
}
dev.off()

logZ
dd
