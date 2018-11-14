# This uses IPTW to estimate the complier stochastic direct effect and its variance, from the following publication
# Kara E. Rudolph, Oleg Sofrygin, and Mark J. van der Laan. "Complier stochastic direct effects: identification and robust estimation." arXiv preprint arXiv:1810.12452 (2018).
# It takes the following arguments:
# a is the instrument, 0/1. It is assumed to be exogenous, but the code can be modified to make it conditionally random.
# z is the exposure influenced by the instrument, 0/1. It is a function of a and w
# m is the mediator, 0/1. It is a function of z, w.
# y is the outcome, 0/1, but the code can be modified for any outcome type. It is a function of z, w, m. 
# w is a matrix of covariates
# svywt is a vector of weights to be applied to the data. 
# zmodel is the parametric model for z.
# mmodel is the parametric model for m.
# ymodel is the parametric model for y. 
# qmodel is the parametric model for q.
# gm is the user-specified stochastic intervention on M, conditional on a=0 and w
# za, za1, and za0 are optional arguments that can be included if the user estimates these as part of the stochastic intervention. Otherwise, they are estimated within the function
# uses the constrained regression function if za, za1, and za0 are null

source("constrainedlogis.R")

mediptw<-function(a, z, m, y, w, svywt, zmodel, mmodel, ymodel, qmodel, gm, za=NULL, za1=NULL, za0=NULL){

datw<-w

# estimate p(m | w, z)
mz<-predict(glm(formula=mmodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m))), newdata=data.frame(cbind(datw, z=z)), type="response")
mz0<-predict(glm(formula=mmodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m))), newdata=data.frame(cbind(datw, z=0)), type="response")
mz1<-predict(glm(formula=mmodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m))), newdata=data.frame(cbind(datw, z=1)), type="response")

# estimate p(z | w, a)
if(is.null(za) | is.null(za1) | is.null(za0)){
 zfit<-mle.logreg.constrained(formula(zmodel), data.frame(cbind(datw, a=a, z=z)))

  za0<-predictClogis(cbind(rep(0,nrow(data.frame(datw))), datw), zfit$beta)
  za1<-predictClogis(cbind(rep(1,nrow(data.frame(datw))), datw), zfit$beta)
  za<-predictClogis(data.frame(cbind(a=a, datw)), zfit$beta)
}
else {
  za<-za
  za1<-za1
  za0<-za0
}

pza1<-ifelse(z==1, za1, 1-za1)
pza0<-ifelse(z==1, za0, 1-za0)

# estimate p(a|w,m,z) using previous estimates. Note that p(a|w,m,z) = p(a|w,z) bc of exclusion restriction 
pa1<-(mean(a)*pza1)/(pza1*mean(a) + pza0*mean(1-a))
pa1z0<-(mean(a)*(1-za1))/((1-za1)*mean(a) + (1-za0)*mean(1-a))
pa1z1<-(mean(a)*za1)/(za1*mean(a) + za0*mean(1-a))

  tmpdat<-data.frame(cbind(datw, a=a))

 #make clever covariate
  psm<-(mz*m) + ((1-mz)*(1-m))
  
  tmpdat$wts<-((m*gm + (1-m)*(1-gm))/psm)* svywt
  #component that can't go into the weights
  tmpdat$cc<- (pa1/mean(a))  - ((1-pa1)/mean(1-a))

  tmpdat$ccz0<- (pa1z0/mean(a))  - ((1-pa1z0)/mean(1-a))
  tmpdat$ccz1<- (pa1z1/mean(a))  - ((1-pa1z1)/mean(1-a))
 
  tmpdat$y<-y
  
  psi1<-sum(tmpdat$y * tmpdat$wts * tmpdat$cc)/sum(svywt) 
  eicpsi1<-(tmpdat$cc*tmpdat$wts * tmpdat$y) - psi1

  #estimate denominator 
  psi2<-sum(z * tmpdat$cc *svywt)/sum(svywt)
  eicpsi2<-(z * tmpdat$cc * svywt) - psi2

  csde<-psi1/psi2
  csdeeic<-(eicpsi1/psi2) - ((psi1*eicpsi2)/(psi2^2))
  varcsde<-var(csdeeic)/nrow(tmpdat)

 return(list("est"=csde, "var"=varcsde))
}
