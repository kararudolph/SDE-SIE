set.seed(2350) 

n<-100000

w0<-rbinom(n, 1, .5)
w1<-rbinom(n,1, .4 + (.2*w0))

probsel<-plogis(-1+ log(4)*w1 + log(4)*w0)
psel<-rbinom(n, 1, probsel)
svywt<-mean(probsel)/probsel

#instrument
a<-rbinom(n, 1, .5)
  
#exposure
z0<-rbinom(n,1,plogis(       - log(2)*w1))
z1<-rbinom(n,1,plogis(log(4) - log(2)*w1))
z<-ifelse(a==1, z1, z0)

#mediator
m0<-rbinom(n, 1, plogis(-log(3)          - log(1.4)*w1))
m1<-rbinom(n, 1, plogis(-log(3) + log(10)- log(1.4)*w1))
m<-ifelse(z==1, m1, m0)

#outcomes
y<-rbinom(n,1, plogis(log(1.2)  + (log(3)*z)  + log(3)*m - log(1.2)*w1 + log(1.2)*w1*z) )

dat<-data.frame(w1=w1, a=a, z=z, m=m, y=y, psel=psel, svywt=svywt, radid_person=seq(1,n,1))
obsdat<-dat[dat$psel==1,]
 
zmodel<-"z ~ a + w1 "
mmodel<-"m ~ z + w1"
ymodel<-"y ~ m + z*w1"
qmodel<-"w1"

#make gm
zfit<-glm(formula=zmodel, family="binomial", data=obsdat)
mfit<-glm(formula=mmodel, family="binomial", data=obsdat)

za0<-predict(zfit, newdata=data.frame(w1=obsdat$w1, a=0), type="response")
za1<-predict(zfit, newdata=data.frame(w1=obsdat$w1, a=1), type="response")

mz1<-predict(mfit, newdata=data.frame(w1=obsdat$w1, z=1), type="response")
mz0<-predict(mfit, newdata=data.frame(w1=obsdat$w1, z=0), type="response")

gm<-(mz1*za0) + (mz0*(1-za0))
gma1<-(mz1*za1) + (mz0*(1-za1))

res<-ivmedtmle(a=obsdat$a, z=obsdat$z, m=obsdat$m, y=obsdat$y, w=data.frame(w1=obsdat$w1), svywt=obsdat$svywt, zmodel=zmodel, mmodel=mmodel, ymodel=ymodel, qmodel=qmodel, gm=gm, gma1=gma1)
