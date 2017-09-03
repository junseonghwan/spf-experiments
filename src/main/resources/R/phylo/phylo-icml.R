rm(list=ls())
dir<-paste("~/Dropbox/Research/repo/spf-experiments/output/icml/", sep="")
folders1<-c("1", "2", "4")
folders2<-c("10", "20", "160")
method_type<-"spf"
N<-length(folders1)
num_rep<-1:4
idx<-4
dd<-matrix(0, ncol=length(folders1), nrow = length(folders2))
hh<-matrix(0, ncol=length(folders1), nrow = length(folders2))
logZ<-matrix(0, ncol=length(folders1), nrow = length(folders2))
#pdf("~/temp/icml.pdf",width=6,height=4,paper='special')
for (i in 1:length(folders1))
{
  folder1<-folders1[i]
  for (j in 1:length(folders2))
  {
    folder2<-folders2[j]
    d1k<-read.csv(paste(dir, folder1, "/", folder2, "/output", idx, "/phylo-", method_type, "-heights.csv", sep=""), header=F)
    truth<-read.csv(paste(dir, folder1, "/", folder2, "/output", idx, "/phylo-", method_type, "-height-truth.csv", sep=""), header=F)
    hh[i,j]<-truth$V1
    plot(density(d1k$V1), col='black', xlim=c(min(c(d1k$V1, truth$V1)), max(c(d1k$V1, truth$V1))), main=paste(idx, i))
    abline(v=truth$V1, col='red')
    q1k<-quantile(d1k$V1, probs = c(0.025, 0.975))
    abline(v=q1k, col='black')
    #print(truth$V1 >= q1k[1] & truth$V1 <= q1k[2])
    #print(mean(d1k$V1 < truth$V1))
    
    trueDistMatrix<-read.csv(paste(dir, folder1, "/", folder2, "/output", idx, "/phylo-pairwise-dist-truth-spf.csv", sep=""), header=T)
    distMatrix<-read.csv(paste(dir, folder1, "/", folder2, "/output", idx, "/phylo-pairwise-dist-spf.csv", sep=""), header=T)
    plot(hclust(dist(trueDistMatrix/2)), main="true tree")
    plot(hclust(dist(distMatrix/2)), main="mean tree")
    diff<-trueDistMatrix - distMatrix
    dd[i,j]<-sum(diff[upper.tri(diff)]^2)
    
    logZ[i,j]<-read.csv(paste(dir, folder1, "/", folder2, "/output", idx, "/phylo-spf-logZ.csv", sep=""), header=F)$V1[3]
  }
}
#dev.off()

logZ
dd
