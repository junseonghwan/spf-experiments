library(ggplot2)

numSimul <- 100
T <- 100
coverage_spf<-matrix(0, ncol = T, nrow = numSimul)
coverage_smc<-matrix(0, ncol = T, nrow = numSimul)
for (i in 1:numSimul)
{
  print(i)
  dir_spf<-paste("Google Drive/Research/repo/spf-experiments/output/ricker-spf/simul", i, "/", sep="")
  x_spf<-read.csv(paste(dir_spf, "ricker-spf-latent.csv", sep=""), header=F)

  dir_smc<-paste("Google Drive/Research/repo/spf-experiments/output/ricker-smc/simul", i, "/", sep="")
  x_smc<-read.csv(paste(dir_smc, "ricker-smc-latent.csv", sep=""), header=F)

  int_spf<-matrix(0, ncol=2, nrow = T+1)
  int_smc<-matrix(0, ncol=2, nrow = T+1)

  for (t in 0:T)
  {
    pop<-read.csv(paste(dir_spf, "particles/population", t, ".csv", sep=""), header=F)$V1
    int_spf[t+1,]<-c(quantile(pop, 0.025), quantile(pop, 0.975))
    pop<-read.csv(paste(dir_smc, "particles/population", t, ".csv", sep=""), header=F)$V1
    int_smc[t+1,]<-c(quantile(pop, 0.025), quantile(pop, 0.975))
  }
  coverage_spf[i,] <- (x_spf >= int_spf[,1] & x_spf <= int_spf[,2])[1:T]
  coverage_smc[i,] <- (x_smc >= int_smc[,1] & x_smc <= int_smc[,2])[1:T]
}

cbind(colMeans(coverage_spf), colMeans(coverage_smc))
cbind(rowMeans(coverage_spf), rowMeans(coverage_smc))

mean(rowMeans(coverage_spf) <= rowMeans(coverage_smc))
mean(colMeans(coverage_spf) <= colMeans(coverage_smc))

# the coverage probability is worse than SMC with larger number of particles
