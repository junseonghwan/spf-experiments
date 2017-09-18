library(ggplot2)

xx<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-sis/simul1/ricker-latent.csv", header=F)$V1
yy<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-sis/simul1/ricker-data.csv", header=F)$V1

# plot of the data
T<-length(yy)
df1<-data.frame(r=1:T, val=yy, type="Observation")
df2<-data.frame(r=1:T, val=xx, type="Latent")
df<-rbind(df1, df2)
p <- ggplot(df, aes(r, val, col=factor(type))) + geom_line() + theme_bw() + theme(legend.title=element_blank())
p <- p + ylab("") + xlab("Time")
p <- ggplot(df1, aes(r, val)) + geom_line() + theme_bw() + theme(legend.position = "none")
p <- p + ylab("Observation") + xlab("Time")
ggsave(filename = "Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-data.pdf", p, width = 4, height = 2)
p2 <- ggplot(df2, aes(r, val)) + geom_line() + theme_bw() + theme(legend.position = "none")
p2 <- p2 + ylab("Latent Variable") + xlab("Time")
p2
ggsave(filename = "Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-latent.pdf", p2, width = 4, height = 2)

# produce the error from filtering mean and the truth as well as ess 
T<-51
err_sis<-rep(0, T)
xhat_sis<-rep(0, T)
err_smc<-rep(0, T)
xhat_smc<-rep(0, T)
for (tt in 0:(T-1))
{
  particles<-read.csv(paste("~/Dropbox/Research/repo/spf-experiments/output/ricker-sis/simul1/population/particles", tt, ".csv", sep=""), header=F)$V1
  weights<-read.csv(paste("~/Dropbox/Research/repo/spf-experiments/output/ricker-sis/simul1/population/weights", tt, ".csv", sep=""), header=F)$V1
  xhat_sis[tt+1]<-sum(weights*particles)
  err_sis[tt+1]<-abs(xx[tt+1] - xhat_sis[tt+1])
  
  particles<-read.csv(paste("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/population/particles", tt, ".csv", sep=""), header=F)$V1
  weights<-read.csv(paste("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/population/weights", tt, ".csv", sep=""), header=F)$V1
  xhat_smc[tt+1]<-sum(weights*particles)
  err_smc[tt+1]<-abs(xx[tt+1] - xhat_smc[tt+1])
}
ess_sis<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-sis/simul1/ricker-ess.csv", header=F)
sis_time<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-sis/simul1/ricker-timing.csv", header=F)
df_sis<-data.frame(Time=0:(T-1), Error=err_sis, ESS=ess_sis$V1, Type="SIS")
p1<-ggplot(df_sis, aes(Time, Error)) + geom_line() + ylab("Absolute Error") + theme_bw()
p2<-ggplot(df_sis, aes(Time, ESS)) + geom_line() + ylab("Effective Sample Size") + theme_bw()
ggsave(p1, filename = "~/Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-sis-error.pdf", width=4, height=2)
ggsave(p2, filename = "~/Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-sis-ess.pdf", width=4, height=2)

ess_smc<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/ricker-ess.csv", header=F)
smc_time<-read.csv("~/Dropbox/Research/repo/spf-experiments/output/ricker-smc/simul1/ricker-timing.csv", header=F)
df_smc<-data.frame(Time=0:(T-1), Error=err_smc, ESS=ess_smc$V1, Type="BPF")
p1<-ggplot(df_smc, aes(Time, Error)) + geom_line() + ylab("Absolute Error") + theme_bw()
p2<-ggplot(df_smc, aes(Time, ESS)) + geom_line() + ylab("Effective Sample Size") + theme_bw()
ggsave(p1, filename = "~/Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-smc-error.pdf", width=4, height=2)
ggsave(p2, filename = "~/Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-smc-ess.pdf", width=4, height=2)

# combine sis and smc plots together
df<-rbind(df_sis, df_smc)
p1<-ggplot(df, aes(Time, Error, col=Type)) + geom_line() + ylab("Absolute Error") + theme_bw() + theme(legend.title=element_blank())
p2<-ggplot(df, aes(Time, ESS, col=factor(Type))) + geom_line() + ylab("Effective Sample Size") + theme_bw() + theme(legend.title=element_blank())
p1
p2
ggsave(p1, filename = "~/Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-error-sis-smc.pdf", width=6, height=2)
ggsave(p2, filename = "~/Dropbox/Research/repo/phd-thesis/contents/smc/figs/ricker-ess-sis-smc.pdf", width=6, height=2)
