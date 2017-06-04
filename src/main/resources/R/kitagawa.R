dsmc<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawaSMC.csv", header=F)
dspf<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawaSPF.csv", header=F)
truth<-read.csv("Google Drive/Research/repo/spf-experiments/output/kitagawaTruth.csv", header=T)
R <- dim(truth)[1]

#plot(density(dsmc$V1), xlim = c(-20, 20))
plot(density(dsmc$V1), ylim = c(0, 0.2))
lines(density(dspf$V1), col='blue', lwd=2, lty=2)
points(truth$x[R], 0, col='red', pch=19)
