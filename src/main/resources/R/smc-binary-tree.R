library(ggplot2)

smcI<-read.csv("Dropbox/Research/repo/spf-experiments/output/testBinaryTree/estimates-incorrect.csv", header=T)

df1<-data.frame(N=smcI$nParticles, mu=smcI$mu_balanced, std=smcI$sd_balanced, type="A balanced Tree")
df2<-data.frame(N=smcI$nParticles, mu=smcI$mu_unbalanced, std=smcI$sd_unbalanced, type="An unbalanced Tree")
df<-rbind(df1, df2)

p<-ggplot(df, aes(N, mu, col=type)) + geom_line() + geom_errorbar(data=df, aes(ymin = mu - std, ymax = mu + std))
p<-p+theme_bw()+xlab("Num. Particles")+ylab("Estimate")+theme(legend.title = element_blank())
p<-p+theme(legend.justification=c(1,1), legend.position=c(1,1))
p<-p+theme(legend.background = element_rect(fill="white", linetype = "dashed"), legend.box.background = element_rect(size=1))
p
ggsave(p, filename = "Dropbox/Research/repo/phd-thesis/contents/smc/figs/smc-binary-tree-incorrect.pdf", width=5.5, height = 5)

rm(list=ls())
smcC<-read.csv("Dropbox/Research/repo/spf-experiments/output/testBinaryTree/estimates-correct.csv", header=T)

df1<-data.frame(N=smcC$nParticles, mu=smcC$mu_balanced, std=smcC$sd_balanced, type="A balanced Tree")
df2<-data.frame(N=smcC$nParticles, mu=smcC$mu_unbalanced, std=smcC$sd_unbalanced, type="An unbalanced Tree")
df<-rbind(df1, df2)

p<-ggplot(df, aes(N, mu, col=type)) + geom_line() + geom_errorbar(data=df, aes(ymin = mu - std, ymax = mu + std))
p<-p+theme_bw()+xlab("Num. Particles")+ylab("Estimate")+theme(legend.title = element_blank())
p<-p+theme(legend.justification=c(1,1), legend.position=c(1,1))
p<-p+theme(legend.background = element_rect(fill="white", linetype = "dashed"), legend.box.background = element_rect(size=1))
p
ggsave(p, filename = "Dropbox/Research/repo/phd-thesis/contents/smc/figs/smc-binary-tree-correct.pdf", width=5.5, height = 5)
