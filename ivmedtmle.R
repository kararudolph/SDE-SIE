#this function uses TMLE to estimate the stochastic direct and indirect effect 
#for observed data O=(W,A,Z,M,Y), where
#W are covariates
#A is an instrumental variable (randomly assigned)
#Z is an intermediate variable that is affected by instrument A
#M is a mediator variable and a function of W, Z (but not A, adhering to exclusion restriction)
#Y is an outcome variable and a function of W, Z, M (but not A, adhering to exclusion restriction)

ivmedtmle<-function(a, z, m, y, w, svywt, zmodel, mmodel, ymodel, qmodel, gm, gma1){
  set.seed(34059)
  datw<-w
  tmpdat<-data.frame(cbind(datw, a=a))

  #get inital fit Q_Y
  tmpdat$qyinit<-cbind(predict(glm(formula=ymodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m, y=y))), newdata=data.frame(cbind(datw, z=z, m=m)), type="response"), 
    predict(glm(formula=ymodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m, y=y))), newdata=data.frame(cbind(datw, z=z, m=0)), type="response"),
    predict(glm(formula=ymodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m, y=y))), newdata=data.frame(cbind(datw, z=z, m=1)), type="response"))
  
  #estimate weights for targeting
  psa1<-I(a==1)/mean(a)
  psa0<-I(a==0)/mean(1-a)
  mz<-predict(glm(formula=mmodel, family="binomial", data=data.frame(cbind(datw, z=z, m=m))), newdata=data.frame(cbind(datw, z=z)), type="response")
  psm<-(mz*m) + ((1-mz)*(1-m))
  
  tmpdat$ha1gma1<-((m*gma1 + (1-m)*(1-gma1))/psm) * psa1 * svywt
  tmpdat$ha1gma0<-((m*gm + (1-m)*(1-gm))/psm) * psa1 * svywt
  tmpdat$ha0gma0<-((m*gm + (1-m)*(1-gm))/psm) * psa0 * svywt
  
  tmpdat$y<-y

  #target Q_Y
  #for E(Y_{1,gmastar})
  epsilonma1g0<-coef(glm(y ~  1 , weights=tmpdat$ha1gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a1g0<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilonma1g0)
  tmpdat$qyupm1a1g0<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilonma1g0)
  #for E(Y_{1,gma})
  epsilonma1g1<-coef(glm(y ~  1 , weights=tmpdat$ha1gma1, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a1g1<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilonma1g1)
  tmpdat$qyupm1a1g1<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilonma1g1)
  #for E(Y_{0,gmastar})
  epsilonma0g0<-coef(glm(y ~  1 , weights=tmpdat$ha0gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a0g0<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilonma0g0)
  tmpdat$qyupm1a0g0<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilonma0g0)

  #estimate Q_M
  tmpdat$Qma1g0<-tmpdat$qyupm0a1g0*(1-gm) + tmpdat$qyupm1a1g0*gm
  tmpdat$Qma1g1<-tmpdat$qyupm0a1g1*(1-gma1) + tmpdat$qyupm1a1g1*gma1
  tmpdat$Qma0g0<-tmpdat$qyupm0a0g0*(1-gm) + tmpdat$qyupm1a0g0*gm

  #estimate Q_Z
  Qzfita1g0<-glm(formula=paste("Qma1g0", qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
  Qzfita1g1<-glm(formula=paste("Qma1g1", qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
  Qzfita0g0<-glm(formula=paste("Qma0g0", qmodel, sep="~"), data=tmpdat[tmpdat$a==0,], family="quasibinomial")
  
  Qza1g0<-predict(Qzfita1g0, type="response", newdata=tmpdat)
  Qza1g1<-predict(Qzfita1g1, type="response", newdata=tmpdat)
  Qza0g0<-predict(Qzfita0g0, type="response", newdata=tmpdat)

  #update Q_Z 
  #Note: only need to do the update step if A is nonrandom
  epsilonza1g0<-coef(glm(Qma1g0~ 1 , weights=psa1*svywt, offset=qlogis(Qza1g0), family="quasibinomial", data=tmpdat))
  epsilonza1g1<-coef(glm(Qma1g1~ 1 , weights=psa1*svywt, offset=qlogis(Qza1g1), family="quasibinomial", data=tmpdat))
  epsilonza0g0<-coef(glm(Qma0g0~ 1 , weights=psa0*svywt, offset=qlogis(Qza0g0), family="quasibinomial", data=tmpdat))

  Qzupa1g0<-plogis(qlogis(Qza1g0) + epsilonza1g0)
  Qzupa1g1<-plogis(qlogis(Qza1g1) + epsilonza1g1)
  Qzupa0g0<-plogis(qlogis(Qza0g0) + epsilonza0g0)
  
  #estimate psi
  tmlea1m0<-sum(Qzupa1g0*svywt)/sum(svywt)
  tmlea1m1<-sum(Qzupa1g1*svywt)/sum(svywt)
  tmlea0m0<-sum(Qzupa0g0*svywt)/sum(svywt)

  #EIC
  tmpdat$qyupa1g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilonma1g0)
  tmpdat$qyupa1g1<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilonma1g1)
  tmpdat$qyupa0g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilonma0g0)

  eic1a1g0<-tmpdat$ha1gma0 * (tmpdat$y - tmpdat$qyupa1g0)
  eic2a1g0<-psa1*svywt*(tmpdat$Qma1g0- Qzupa1g0)
  eic3a1g0<-Qzupa1g0-tmlea1m0
  eica1g0<-eic1a1g0 + eic2a1g0 + eic3a1g0

  eic1a1g1<-tmpdat$ha1gma1 * (tmpdat$y - tmpdat$qyupa1g1)
  eic2a1g1<-psa1*svywt*(tmpdat$Qma1g1- Qzupa1g1)
  eic3a1g1<-Qzupa1g1-tmlea1m1
  eica1g1<-eic1a1g1 + eic2a1g1 + eic3a1g1

  eic1a0g0<-tmpdat$ha0gma0 * (tmpdat$y - tmpdat$qyupa0g0)
  eic2a0g0<-psa0*svywt*(tmpdat$Qma0g0- Qzupa0g0)
  eic3a0g0<-Qzupa0g0-tmlea0m0
  eica0g0<-eic1a0g0 + eic2a0g0 + eic3a0g0

  #estimands
  nde<-tmlea1m0-tmlea0m0
  ndeeic<-eica1g0 - eica0g0
  vareic<-var(ndeeic)/nrow(tmpdat)

  nie<-tmlea1m1-tmlea1m0
  nieeic<-eica1g1 - eica1g0
  varnieeic<-var(nieeic)/nrow(tmpdat)

  return(list("nde"=nde, "ndevar"=vareic, "nie"=nie, "nievar"=varnieeic))
}
