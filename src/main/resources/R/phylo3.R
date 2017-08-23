data_dir <- "~/Google Drive/Research/repo/spf-experiments/output/"
truth_smc<-read.csv(paste(data_dir, "/phylo-smc-height-truth.csv", sep=""), header=F)$V1
samples_smc<-read.csv(paste(data_dir, "/phylo-smc-heights.csv", sep=""), header=F)$V1
weights_smc<-read.csv(paste(data_dir, "/phylo-smc-weights.csv", sep=""), header=F)$V1

truth_spf<-read.csv(paste(data_dir, "/phylo-spf-height-truth.csv", sep=""), header=F)$V1
samples_spf<-read.csv(paste(data_dir, "/phylo-spf-heights.csv", sep=""), header=F)$V1

plot(density(x=samples_smc, weights = weights_smc), col='red')
abline(v=truth_smc, col='green')
lines(density(x=samples_spf), lty=2, col='blue')
abline(v=truth_spf, col='black')

# 95% CI
quantile(samples_smc, c(0.025, 0.975))
quantile(samples_spf, c(0.025, 0.975))
truth_spf
