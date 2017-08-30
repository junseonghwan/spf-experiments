library(ggplot2)
library(reshape2)
library(pomp)
library(coda)

pompExample(ricker)
ricker
plot(ricker)
x<-simulate(ricker)
as.data.frame(x)
plot(x)
y<-as.data.frame(ricker)
head(y)
head(simulate(ricker,as.data.frame=TRUE))
x <- simulate(ricker,nsim=10)
class(x)
sapply(x,class)

x <- simulate(ricker,nsim=10,as.data.frame=TRUE)
head(x)
str(x)

x <- simulate(ricker,nsim=9,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,aes(x=time,y=y,group=sim,color=(sim=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)

y <- trajectory(ricker)
dim(y)
dimnames(y)
plot(time(ricker),y["N",1,],type="l")
coef(ricker)

? pfilter
pp<-coef(ricker)
pf <- pfilter(object=x, params = pp, Np=1000, pred.mean = TRUE, pred.var = TRUE, filter.mean = TRUE, filter.traj = TRUE, save.states = TRUE, save.params=TRUE)
pf@pred.mean
pf@filter.mean
df<-as.data.frame(x)
length(pf@saved.states)
for (tt in 1:T)
{
  plot(density(pf@saved.states[[tt]]))
  abline(v=df$N[tt], col='red')
}
mean((pf@pred.mean["N",] - df$N)^2)
mean((pf@filter.mean["N",] - df$N)^2)
plot(pf@filter.mean["N",], type='l')
lines(df$N, col='red')
df$y
logLik(pf)
pf@loglik

pf <- pfilter(pf)
logLik(pf)

# run pmcmc
x
chain<-pmcmc(x, Nmcmc = 1000, Np = 1000, verbose=TRUE, 
             start=c("r"=44,"sigma"=0.2,"phi"=9,"N.0"=7,"e.0"=0), 
             proposal=mvn.rw.adaptive(rw.sd=c(r=0.01,sigma=0.01,phi=0.01, N.0=0.01), scale.start=200,shape.start=100))
plot(chain)
chain2<-pmcmc(chain)
trace <- window(conv.rec(chain,c("r","sigma","phi","N.0")),start=100)
rejectionRate(trace)
effectiveSize(trace)
autocorr.diag(trace)
summary(trace)
plot(trace)
heidel.diag(trace)
geweke.diag(trace)
chain@conv.rec
#bsflu_data <- read.table("https://kingaa.github.io/sbied/stochsim/bsflu_data.txt")
