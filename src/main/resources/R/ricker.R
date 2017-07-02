library(ggplot2)
x<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-smc-latent.csv", header=F)$V1
xx<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-spf-latent.csv", header=F)$V1

y<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-smc-data.csv", header=F)$V1
yy<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-spf-data.csv", header=F)$V1

prior_samples<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-prior-samples.csv", header=F)
dim(prior_samples)
plot(x=0,y=0, type='n', xlim = c(0, 11), ylim=c(0, max(prior_samples)))
for (i in 1:100)
{
  lines(0:10, prior_samples[i,], type='l')
}

neff<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-smc-ess.csv", header=F)
nimp<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-spf-nimplicit.csv", header=F)
smc_time<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-smc-timing.csv", header=F)
spf_time<-read.csv("~/Google Drive/Research/repo/spf-experiments/output/ricker-spf-timing.csv", header=F)
sum(smc_time$V1)
sum(spf_time$V1)

T<-1001
err_spf<-rep(0, T)
err_smc<-rep(0, T)
var_spf<-rep(0, T)
var_smc<-rep(0, T)
mu_spf<-rep(0, T)
mu_smc<-rep(0, T)
interval_smc<-matrix(0, ncol = 2, nrow = T)
interval_spf<-matrix(0, ncol = 2, nrow = T)
for (t in 0:(T-1))
{
  dd<-read.csv(paste("~/Google Drive/Research/repo/spf-experiments/output/ricker-spf/particles", t, ".csv", sep=""), header=F)
  #d<-read.csv(paste("~/Google Drive/Research/repo/spf-experiments/output/ricker-smc/particles", t, ".csv", sep=""), header=F)
  
  #interval_smc[t+1,]<-quantile(d$V1, probs = c(0.025, 0.975))
  interval_spf[t+1,]<-quantile(dd$V1, probs = c(0.025, 0.975))
  
  err_spf[t+1]<-mean((dd$V1 - x[t+1])^2)
  #err_smc[t+1]<-mean((d$V1 - x[t+1])^2)
  
  mu_spf[t+1]<-mean(dd$V1)
  #mu_smc[t+1]<-mean(dd$V1)
  
  var_spf[t+1]<-var(dd$V1)
  #var_smc[t+1]<-var(d$V1)
  
  # p <- ggplot(dd, aes(V1)) + geom_density() + theme_bw()
  # p <- p + geom_point(x=x[t+1], y=0, col='red')
  # ggsave(paste("Google Drive/Research/repo/spf-experiments/output/ricker-plots/density-spf-", t, ".pdf", sep=""), p)
  # 
  # p <- ggplot(d, aes(V1)) + geom_density() + theme_bw()
  # p <- p + geom_point(x=x[t+1], y=0, col='red')
  # ggsave(paste("Google Drive/Research/repo/spf-experiments/output/ricker-plots/density-smc-", t, ".pdf", sep=""), p)
}
df1<-data.frame("t"=1:T, "err"=err_spf, type="spf")
df2<-data.frame("t"=1:T, "err"=err_smc, type="smc")
df<-rbind(df1, df2)
ggplot(df, aes(t, err ,col=type)) + geom_point(position="jitter")

plot(1:T, neff$V1)
plot(1:T, nimp$V1)

plot(1:T, var_smc, col='red', pch=19)
points(1:T, var_spf, col='blue', pch=19)
boxplot(var_smc)
boxplot(var_spf)

plot(1:T, y, col='blue', pch=19)
points(1:T, 10*x, col='red', pch=19)
points(1:T, 10*mu_spf, col='black', pch=19)

plot(1:T, y, col='blue', pch=19)
points(1:T, 10*x, col='red', pch=19)
points(1:T, 10*mu_smc, col='orange', pch=17)

# check the coverage of the 95% confidence region
mean(x >= interval_smc[,1] & x <= interval_smc[,2])
mean(x >= interval_spf[,1] & x <= interval_spf[,2])

# length of the CI
l_smc <- interval_smc[,2] - interval_smc[,1]
l_spf <- interval_spf[,2] - interval_spf[,1]
par(mfrow=c(1,2))
boxplot(l_smc)
boxplot(l_spf)

mean(l_smc < l_spf)

# since the data is generated from the prior, posterior should converge to prior as the data size increases
# so, check how many of the samples from the prior are in the CI built by the particles
dim(interval_spf)
dim(prior_samples)
for (i in 1:dim(prior_samples)[1])
{
  mean(t(prior_samples[i,]) >= interval_spf[,1] & t(prior_samples[i,]) <= interval_spf[,2])
}

plot(x=0,y=0, type='n', xlim = c(0, T-1), ylim=c(0, max(prior_samples)))
for (i in 1:100)
{
  points(0:(T-1), prior_samples[i,])
}

### PMMH results

d<-read.csv("Google Drive/Research/repo/spf-experiments/output/2017-06-29-06-39-15-Oy7C0d7r.exec/spf.pmcmc.PMCMCDefaultOutputProcessor.csv", header=F)
pairs(d)
truth<-c(10, 44.7, 0.3)
plot(d$V1, d$V2)
points(truth[1], truth[2], col='red', pch=19)
plot(d$V1, d$V3)
points(truth[1], truth[3], col='red', pch=19)
plot(d$V2, d$V3)
points(truth[2], truth[3], col='red', pch=19)
plot(d$V1, type='l')
abline(h=truth[1], col='red')
plot(d$V2, type='l')
abline(h=truth[2], col='red')
plot(d$V3, type='l')
abline(h=truth[3], col='red')
