rm(list=ls())
method_type<-"spf"
dir1k<-paste("~/Dropbox/Research/repo/spf-experiments/output/spf20taxa/", sep="")
numSimul<-1
sum1k <- rep(0, numSimul)
quan1k<-rep(0, numSimul)
hh<-rep(0, numSimul)
dd<-rep(0, numSimul)
pdf("~/temp/file1.pdf",width=6,height=4,paper='special')
for (i in 1:numSimul)
{
  d1k<-read.csv(paste(dir1k, "output", i, "/phylo-", method_type, "-heights.csv", sep=""), header=F)
  truth<-read.csv(paste(dir1k, "output", i, "/phylo-", method_type, "-height-truth.csv", sep=""), header=F)
  hh[i]<-truth$V1
  plot(density(d1k$V1), col='black', xlim=c(min(c(d1k$V1, truth$V1)), max(c(d1k$V1, truth$V1))), main = paste("Simulation", i))
  abline(v=truth$V1, col='red')
  q1k<-quantile(d1k$V1, probs = c(0.025, 0.975))
  abline(v=q1k, col='black')
  sum1k[i] <- (truth$V1 >= q1k[1] & truth$V1 <= q1k[2])
  quan1k[i]<-mean(d1k$V1 < truth$V1)
  
  trueDistMatrix<-read.csv(paste(dir1k, "output", i, "/phylo-pairwise-dist-truth-spf.csv", sep=""), header=T)
  distMatrix<-read.csv(paste(dir1k, "output", i, "/phylo-pairwise-dist-spf.csv", sep=""), header=T)
  plot(hclust(dist(trueDistMatrix)), main="true tree")
  plot(hclust(dist(distMatrix)), main="mean tree")
  diff<-trueDistMatrix - distMatrix
  dd[i]<-mean(diff[upper.tri(trueDistMatrix)]^2)
  plot(density(diff[upper.tri(trueDistMatrix)]))
}
dev.off()
sum(sum1k)/numSimul
which(sum1k == 0)

hist(quan1k, breaks = 10)
z<-qnorm(quan1k)
which(z == Inf)
which(z == -Inf)
zz<-z[abs(z) != Inf]
hist(zz)
plot(density(zz))
c(mean(zz), sd(zz))
hist(zz^2, breaks=20)
plot(1:length(zz), abs(zz))
qq<-sum(zz^2)
qq
pchisq(qq, df = length(zz), lower.tail = TRUE)
xx<-rchisq(100000, df = length(zz))
plot(density(xx))
abline(v=qq, col='red')
mean(qq >= xx)
