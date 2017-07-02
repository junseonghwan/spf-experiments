dsmc<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawaSMC.csv", header=F)
dspf<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawaSPF.csv", header=F)
truth<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawaTruth.csv", header=T)
T <- dim(truth)[1]

#plot(density(dsmc$V1), xlim = c(-20, 20))
plot(density(dsmc$V1), ylim = c(0, 0.2))
lines(density(dspf$V1), col='blue', lwd=2, lty=2)
points(truth$x[R], 0, col='red', pch=19)

err_spf<-rep(0, T)
err_smc<-rep(0, T)
var_spf<-rep(0, T)
var_smc<-rep(0, T)
mu_spf<-rep(0, T)
mu_smc<-rep(0, T)

for (t in 1:T)
{
  smc<-read.csv(paste("Google Drive/Research/repo/spf-experiments/output/kitagawa-smc/particles", (t-1), ".csv", sep=""), header=F)
  spf<-read.csv(paste("Google Drive/Research/repo/spf-experiments/output/kitagawa-spf/particles", (t-1), ".csv", sep=""), header=F)
  
  err_spf[t]<-mean((spf$V1 - truth$x[t])^2)
  err_smc[t]<-mean((smc$V1 - truth$x[t])^2)

  var_spf[t]<-var(spf$V1)
  var_smc[t]<-var(smc$V1)

  mu_spf[t]<-mean(spf$V1)
  mu_smc[t]<-mean(smc$V1)
  
  #plot(density(smc$V1))
  #lines(density(spf$V1), lty=2, col='blue')
  #points(x=truth$x[t], y=0, col='red')
}

