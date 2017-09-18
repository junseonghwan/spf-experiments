library(reshape2)
library(ggplot2)
library(dplyr)

plot_w<-5
plot_h<-3

f<-function(nTaxa, K, method_type, numrep=10)
{
  base_dir<-"~/Dropbox/Research/repo/spf-experiments/output/phylo-timing/"
  dir<-paste(base_dir, "nTaxa", nTaxa, "/", method_type, "/output", sep="")
  timing<-matrix(0, nrow = numrep, ncol = nTaxa-1)
  rawESS<-matrix(0, nrow = numrep, ncol = nTaxa-1)
  essPerSec<-matrix(0, nrow = numrep, ncol = nTaxa-1)
  for (output_idx in 1:numrep)
  {
    dir2<-paste(dir, output_idx, sep="")
    d<-read.csv(paste(dir2, "/pairwise-dist-", method_type, ".csv", sep=""), header=T)
    d_truth<-read.csv(paste(dir2, "/pairwise-dist-truth.csv", sep=""), header=T)
    
    h<-read.csv(paste(dir2, "/heights-", method_type, ".csv", sep=""), header=F)$V1
    h_truth<-read.csv(paste(dir2, "/height-truth.csv", sep=""), header=F)$V1
    w<-rep(1/K, K)
    if (method_type == "smc")
    {
      w<-read.csv(paste(dir2, "/weights-", method_type, ".csv", sep=""), header=F)$V1
    }
    
    #print(sum((d - d_truth)^2)/2)
    #print(sum(h*w) - h_truth)
    #plot(density(h, weights = w))
    #abline(v=h_truth)
    
    timing[output_idx,]<-read.csv(paste(dir2, "/timing-results-", method_type, ".csv", sep=""), header=F)$V1
    rawESS[output_idx,]<-read.csv(paste(dir2, "/ess-", method_type, ".csv", sep=""), header=F)$V1
    essPerSec[output_idx,]<-rawESS[output_idx,]/timing[output_idx,]
  }
  ret<-list("rawESS"=rawESS, "timing"=timing, "essPerSec"=essPerSec)
  return(ret)
}

# summarize the data frame by iteration
g<-function(dd)
{
  mu<-dd %>% group_by(iter) %>% summarise(value = mean(value))
  sde<-dd %>% group_by(iter) %>% summarise(value = sd(value))
  dd2<-data.frame("iter"=mu$iter, "avg"=mu$value, "std_err"=sde$value)
  return(dd2)
}

nTaxa<-16
K<-10000
xspf<-f(nTaxa, K, "spf")
xsmc<-f(nTaxa, K, "smc")
ratio<-xspf$essPerSec/xsmc$essPerSec

# plot the raw ESS at each iteration
dd1<-g(melt(xspf$rawESS, varnames = c("rep", "iter")))
dd1$type<-"SPF"
dd2<-g(melt(xsmc$rawESS, varnames = c("rep", "iter")))
dd2$type<-"SMC"
dd<-rbind(dd1, dd2)
p<-ggplot(dd, aes(iter, avg, col=type)) + geom_point()
p<-p+geom_errorbar(aes(ymin=avg-std_err, ymax=avg+std_err))
p<-p+theme_bw()+ylab("ESS")+xlab("SMC Iteration")
p<-p+theme(legend.justification=c(0,1), legend.position=c(0,1), legend.title = element_blank(), legend.background = element_blank())
p
ggsave(plot = p, filename = paste("Dropbox/Research/repo/phd-thesis/contents/phylogenetics/figures/phylo-ess-", nTaxa, ".pdf", sep=""), width = plot_w, height = plot_h)

# compute the average times
dd1<-g(melt(xspf$timing, varnames = c("rep", "iter")))
dd1$type<-"SPF"
dd2<-g(melt(xsmc$timing, varnames = c("rep", "iter")))
dd2$type<-"SMC"
dd<-rbind(dd1, dd2)
dd$avg<-dd$avg*1000
dd$std_err<-dd$std_err*1000
p<-ggplot(dd, aes(iter, avg, col=type)) + geom_point()
p<-p+geom_errorbar(aes(ymin=avg-std_err, ymax=avg+std_err))
p<-p+theme_bw()+ylab("Time (in milliseconds)")+xlab("SMC Iteration")
p<-p+theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title = element_blank(), legend.background = element_blank())
p
ggsave(plot = p, filename = paste("Dropbox/Research/repo/phd-thesis/contents/phylogenetics/figures/phylo-timing-", nTaxa, ".pdf", sep=""), width = plot_w, height = plot_h)

# plot the ESS/sec
dd1<-g(melt(xspf$essPerSec, varnames = c("rep", "iter")))
dd1$type<-"SPF"
dd2<-g(melt(xsmc$essPerSec, varnames = c("rep", "iter")))
dd2$type<-"SMC"
dd<-rbind(dd1, dd2)
p <- ggplot(dd, aes(iter, avg, col=type)) + geom_line() + geom_point() 
p <- p + theme_bw()
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1), legend.title = element_blank(), legend.background = element_blank())
p <- p + xlab("SMC Iteration") + ylab("ESS/sec")
p
ggsave(plot = p, filename = paste("Dropbox/Research/repo/phd-thesis/contents/phylogenetics/figures/ess-per-sec-", nTaxa, ".pdf", sep=""), width = plot_w, height = plot_h)

# plot the ratio of ESS/sec
dd<-melt(ratio, varnames = c("rep", "iter"))
p <- ggplot(dd, aes(iter, value)) + geom_smooth() + geom_hline(aes(yintercept = 1), col="red")
p <- p + theme_bw() + xlab("SMC Iterations") + ylab("Ratio of ESS/sec")
p
ggsave(plot = p, filename = paste("Dropbox/Research/repo/phd-thesis/contents/phylogenetics/figures/ess-ratio-", nTaxa, ".pdf", sep=""), width = plot_w, height = plot_h)
