d_smc<-read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-pairwise-dist-smc.csv", header=T)
d_spf<-read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-pairwise-dist-spf.csv", header=T)
truth<-read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-pairwise-truth-smc.csv", header=T)

sum(abs(truth - d_smc))/2
sum(abs(truth - d_spf))/2

# read in the ESS
read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-smc-ess.csv", header=F)
read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-spf-ess.csv", header=F)

h_spf <- read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-spf-heights.csv", header=F)
h_smc <- read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-smc-heights.csv", header=F)
w_smc <- read.csv("Google Drive/Research/repo/spf-experiments/output/phylo-smc-weights.csv", header=F)
plot(density(h_spf$V1), lty=2, lwd=2, col='blue', xlim=c(0, max(h_smc)))
lines(density(h_smc$V1, weights = w_smc$V1))
#points(x=1.597091798902487, y = 0, col='red', pch=19)
#abline(v=0.8763986225509377, col='red', lty=2)
abline(v=1.0394309877726353, col='red', lty=2)
