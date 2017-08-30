library(ggplot2)
xx<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/ricker-smc-latent.csv", header=F)$V1
yy<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/ricker-smc-data.csv", header=F)$V1
neff<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/ricker-smc-ess.csv", header=F)
smc_time<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/ricker-smc-timing.csv", header=F)
sum(smc_time$V1)

T<-51
err_smc<-rep(0, T)
var_smc<-rep(0, T)
mu_smc<-rep(0, T)
interval_smc<-matrix(0, ncol = 2, nrow = T)
for (tt in 0:(T-1))
{
  dd<-read.csv(paste("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/particles/population", tt, ".csv", sep=""), header=F)

  interval_smc[tt+1,]<-quantile(dd$V1, probs = c(0.025, 0.975))
  err_smc[tt+1]<-mean((dd$V1 - xx[tt+1])^2)
  mu_smc[tt+1]<-mean(dd$V1)
  var_smc[tt+1]<-var(dd$V1)
  
  # p <- ggplot(dd, aes(V1)) + geom_density() + theme_bw()
  # p <- p + geom_point(x=x[t+1], y=0, col='red')
  # ggsave(paste("Google Drive/Research/repo/spf-experiments/output/ricker-plots/density-spf-", t, ".pdf", sep=""), p)
  # 
  # p <- ggplot(d, aes(V1)) + geom_density() + theme_bw()
  # p <- p + geom_point(x=x[t+1], y=0, col='red')
  # ggsave(paste("Google Drive/Research/repo/spf-experiments/output/ricker-plots/density-smc-", t, ".pdf", sep=""), p)
}
plot(1:T, neff$V1, type='l')
plot(1:T, yy, col='blue', pch=19, type='l')
lines(1:T, xx, col='red', pch=19, type='l')
lines(1:T, mu_smc, col='black', pch=19, type='l')
mean((mu_smc - xx)^2)
