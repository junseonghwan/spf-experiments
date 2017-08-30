dir<-paste("~/Dropbox/Research/repo/spf-experiments/output/icml/", sep="")
folders<-c("10", "80", "160", "320")
method_type<-"spf"
N<-length(folders)
dd<-matrix(0, ncol=N, nrow = 5)
hh<-matrix(0, ncol=N, nrow = 5)
logZ<-matrix(0, ncol=N, nrow = 5)
pdf("~/temp/icml.pdf",width=6,height=4,paper='special')
for (idx in 1:5)
{
  for (i in 1:N)
  {
    folder<-folders[i]
    d1k<-read.csv(paste(dir, folder, "/output", idx, "/phylo-", method_type, "-heights.csv", sep=""), header=F)
    truth<-read.csv(paste(dir, folder, "/output", idx, "/phylo-", method_type, "-height-truth.csv", sep=""), header=F)
    hh[idx,i]<-mean(d1k$V1)-truth$V1
    plot(density(d1k$V1), col='black', xlim=c(min(c(d1k$V1, truth$V1)), max(c(d1k$V1, truth$V1))), main=paste(idx, i))
    abline(v=truth$V1, col='red')
    q1k<-quantile(d1k$V1, probs = c(0.025, 0.975))
    abline(v=q1k, col='black')
    print(truth$V1 >= q1k[1] & truth$V1 <= q1k[2])
    print(mean(d1k$V1 < truth$V1))
    
    trueDistMatrix<-read.csv(paste(dir, folder, "/output", idx, "/phylo-pairwise-dist-truth-spf.csv", sep=""), header=T)
    distMatrix<-read.csv(paste(dir, folder, "/output", idx, "/phylo-pairwise-dist-spf.csv", sep=""), header=T)
    plot(hclust(dist(trueDistMatrix)), main="true tree")
    plot(hclust(dist(distMatrix)), main="mean tree")
    diff<-trueDistMatrix - distMatrix
    dd[idx,i]<-sum(diff[upper.tri(diff)]^2)
    plot(density(diff[upper.tri(trueDistMatrix)]))
    
    logZ[idx,i]<-read.csv(paste(dir, folder, "/output", idx, "/phylo-spf-logZ.csv", sep=""), header=F)$V1
  }
}
dev.off()

logZ
dd
