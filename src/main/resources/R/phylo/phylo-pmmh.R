rm(list=ls())
library(ggplot2)
f<-function(method_type, thinning)
{
  base_dir<-"~/Dropbox/Research/repo/spf-experiments/output/phylo-pmmh/"
  data_file<-paste(method_type, "/", method_type, ".pmcmc.PMCMCDefaultOutputProcessor.csv", sep="")
  x<-read.csv(paste(base_dir, data_file, sep=""), header=F)
  xx<-x[seq(1, dim(x)[1], thinning),]
  return(xx)
}

g<-function(method_type)
{
  base_dir<-"~/Dropbox/Research/repo/spf-experiments/output/phylo-pmmh/"
  data_file<-paste(method_type, "/smcStat.csv", sep="")
  x<-read.csv(paste(base_dir, data_file, sep=""), header=T)
  return(x)
}

plot_w<-5
plot_h<-3

xspf<-f("spf", 30)
xsmc<-f("smc", 30)
df1<-data.frame(val=xspf, type="SPF")
df2<-data.frame(val=xsmc, type="SMC")
df<-rbind(df1, df2)

p <- ggplot(df, aes(val, col=type)) + geom_density() + theme_bw()
p <- p + geom_vline(aes(xintercept = 1.7))
p <- p + xlab("Parameter values") + ylab("Density")
p <- p + theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title = element_blank(), legend.background = element_blank())
p
ggsave(p, filename = "~/Dropbox/Research/repo/phd-thesis/contents/phylogenetics/figures/pmmh-spf.pdf", width = plot_w, height = plot_h)

# read the smc stat
xxspf<-g("spf")
xxsmc<-g("smc")
R<-dim(xxspf)[1]
df1<-data.frame("Iter"=1:R, "val"=xxspf$Avg/xxspf$TimeAvg, type="SPF")
df2<-data.frame("Iter"=1:R, "val"=xxsmc$Avg/xxsmc$TimeAvg, type="SMC")
df<-rbind(df1, df2)
p <- ggplot(df, aes(Iter, val, col=type)) + geom_line() + geom_point()
p <- p + scale_y_log10() + theme_bw()
p <- p + xlab("SMC Iteration") + ylab("ESS/sec (log scale)")
p <- p + theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title = element_blank(), legend.background = element_blank())
p  
ggsave(p, filename = "~/Dropbox/Research/repo/phd-thesis/contents/phylogenetics/figures/pmmh-ess-per-sec.pdf", width= plot_w, height = plot_h)
