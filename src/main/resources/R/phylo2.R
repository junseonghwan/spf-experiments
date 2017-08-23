smc_dir<-"~/Google Drive/Research/repo/spf-experiments/output/phylo-smc/"
spf_dir<-"~/Google Drive/Research/repo/spf-experiments/output/phylo-spf/"
numSimul<-20
x_spf<-rep(0, numSimul)
x_smc<-rep(0, numSimul)
smc_i<-matrix(0, ncol=2, nrow=numSimul)
spf_i<-matrix(0, ncol=2, nrow=numSimul)
for (i in 1:numSimul)
{
  truth_smc<-read.csv(paste(smc_dir, "output", i, "/phylo-smc-height-truth.csv", sep=""), header=F)$V1
  samples_smc<-read.csv(paste(smc_dir, "output", i, "/phylo-smc-heights.csv", sep=""), header=F)$V1
  weights_smc<-read.csv(paste(smc_dir, "output", i, "/phylo-smc-weights.csv", sep=""), header=F)$V1

  truth_spf<-read.csv(paste(spf_dir, "output", i, "/phylo-spf-height-truth.csv", sep=""), header=F)$V1
  samples_spf<-read.csv(paste(spf_dir, "output", i, "/phylo-spf-heights.csv", sep=""), header=F)$V1

  # 95% CI
  smc_i[i,]<-quantile(samples_smc, c(0.025, 0.975))
  spf_i[i,]<-quantile(samples_spf, c(0.025, 0.975))
  x_smc[i]<-truth_smc
  x_spf[i]<-truth_spf
  
  plot(density(samples_spf))
  abline(v=truth_spf, col='red')
}

cbind(x_smc, x_spf)

mean(x_smc >= smc_i[,1] & x_smc <= smc_i[,2])
mean(x_spf >= spf_i[,1] & x_spf <= spf_i[,2])

cbind(spf_i[,1], x_spf, spf_i[,2])

