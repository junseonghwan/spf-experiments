#read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-09-23-51-44-4IQbh373.exec/executionInfo/options.map", sep="\t", header=F)
# SPF
# nConcrete = {100, 200, 500, 1000}
spfNumAccepts<-rep(0, 4)
spfNumAccepts[1]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-09-23-51-44-4IQbh373.exec/summary.txt", sep=":", header=F)[1,2]
spfNumAccepts[2]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-09-23-51-32-1Qvw7Ben.exec/summary.txt", sep=":", header=F)[1,2]
spfNumAccepts[3]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-08-07-42-16-lUTwcl5p.exec/summary.txt", sep=":", header=F)[1,2]
spfNumAccepts[4]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-07-06-28-44-DRI57xJG.exec/summary.txt", sep=":", header=F)[1,2]

# nParticles = {100, 1000, 2000, 5000}
smcNumAccepts<-rep(0, 4)
smcNumAccepts[1]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-06-03-18-17-hxJIZzgy.exec/summary.txt", sep=":", header=F)[1,2]
smcNumAccepts[2]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-06-18-17-08-RdgEzW6O.exec/summary.txt", sep=":", header=F)[1,2]
smcNumAccepts[3]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-06-18-17-41-rnWjSbV1.exec/summary.txt", sep=":", header=F)[1,2]
smcNumAccepts[4]<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-06-23-07-22-t2esy5qZ.exec/summary.txt", sep=":", header=F)[1,2]

plot(spfNumAccepts)
points(smcNumAccepts, col='red')

# read the timing results
df_spf<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-07-06-28-44-DRI57xJG.exec/smcStat.csv", header=T)
df_smc<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-06-23-07-22-t2esy5qZ.exec/smcStat.csv", header=T)
sum(df_spf$TimeAvg)
sum(df_smc$TimeAvg)

df_spf<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-07-06-28-44-DRI57xJG.exec/spf.experiments.kitagawa.KitagawaProcessor.csv", header=F)
df_smc<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawa-pmmh/2017-06-06-23-07-22-t2esy5qZ.exec/smc.experiments.kitagawa.KitagawaProcessor.csv", header=F)
xx_spf<-df_spf[10001:dim(df_spf)[1],]
plot(xx_spf$V1, xx_spf$V2, col='blue', pch=19)

xx_smc<-df_smc[10001:dim(df_smc)[1],]
points(xx_smc$V1, xx_smc$V2, col='red', pch=17)
