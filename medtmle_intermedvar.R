#-----------------------------------------------------------------------------------------------------------------------------
#     
# This is a function to estimate stochastic direct and indirect effects when there are direct effects of the 
#   exposure on the mediator and the outcome, and there is a mediator-outcome confounder affected by prior exposure. 
# This corresponds to the following data-generating structure: 
#     W=f(U_W)
#     A=f(W, U_A)
#     Z=f(W, A, U_Z)
#     M=f(W, A, Z, U_M)
#     Y=f(W, A, Z, M, U_Y)

# The stochastic direct effect is defined as E(Y_{1,g_0}) - E(Y_{0,g_0})
#   and the stochastic indirect effect is defined as E(Y_{1,g_1}) - E(Y_{1,g_0})
#   where Y_{a,g_a*} is the counterfactual value of Y when a=a and the mediator m is drawn from 
#   the distribution g under a*. 
#
# Inputs to the function: 
#     1. obsdat 
#       - should be a data frame that includes an exposure variable (called a with values 0/1), 
#              a mediator-outcome confounder affected by prior exposure (called z with values 0/1), 
#              a mediator variable (called m with values 0/1), and an outcome variable (called y with values 0/1). 
#       - there can be other covariates that can be named anything and have any values.
#     2. amodel 
#       - the parametric model for a that includes the parents of the exposure. For example: "a ~ w1 + w2 + w3"
#     3. zmodel 
#       - the parametric model for z that includes the parents of the mediator-outcome confounder 
#              affected by prior exposure. For example: "z ~ a + w1 + w2 + w3" 
#     4. mmodel 
#       - the parametric model for m that includes the parents of the mediator. 
#     5. ymodel 
#       - the parametric model for y that includes the parents of the outcome. 
#     6. qmodel 
#       - the covariates to marginalize over. For example, "w1 + w2 + w3"
#  
#      
#       
#-----------------------------------------------------------------------------------------------------------------------------

medtmle_intermedvar <- function(obsdat, amodel, zmodel, mmodel, ymodel, qmodel) {
  
  dfa1 <- dfa0 <- dfa1z1 <- dfa1z0 <- dfa0z1 <- dfa0z0 <- dfm1 <- dfm0 <- tmpdat <- obsdat
  
  dfa1$a <- dfa1z1$a <- dfa1z1$z <- dfa1z0$a <- dfa0z1$z <- dfm1$m <- 1 
  dfa0$a <- dfa1z0$z <- dfa0z1$a <- dfa0z0$a <- dfa0z0$z <- dfm0$m <- 0
  

  # calculate gm under a=0 and gm under a=1, then marginalize over z 
  za0 <- predict(glm(formula=zmodel, family="binomial", data=obsdat), newdata=dfa0, type="response")
  za1 <- predict(glm(formula=zmodel, family="binomial", data=obsdat), newdata=dfa1, type="response")
  
  mz1a0 <- predict(glm(formula=mmodel, family="binomial", data=obsdat), newdata=dfa0z1, type="response")
  mz0a0 <- predict(glm(formula=mmodel, family="binomial", data=obsdat), newdata=dfa0z0, type="response")
  
  mz1a1 <- predict(glm(formula=mmodel, family="binomial", data=obsdat), newdata=dfa1z1, type="response")
  mz0a1 <- predict(glm(formula=mmodel, family="binomial", data=obsdat), newdata=dfa1z0, type="response")
  
  # empirical distribution of m under a=0, marginalized over z
  gma0 <- (mz1a0*za0) + (mz0a0*(1-za0))
  # empirical distribution of m under a=1, marginalized over z
  gma1 <- (mz1a1*za1) + (mz0a1*(1-za1))

  # calculate a weights
  psa1 <- I(tmpdat$a==1)/predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response")
  psa0 <- I(tmpdat$a==0)/(1-predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response"))
  
  # calculate m weights
  mazw <-predict(glm(formula=mmodel, family="binomial", data=tmpdat), type="response")
  psm<- I(tmpdat$m==1)*mazw + I(tmpdat$m==0)*(1-mazw)
  
  
 # first generate initial estimate of y in full pop, and then under m=0 and m=1 
  tmpdat$qyinit<-cbind(predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=tmpdat, type="response"), 
                       predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm0, type="response"),
                       predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm1, type="response"))
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_0})
  #-------------------------------------------------------------------------------------------------------------
    # calculate weights when a=1 and gm when a=0 
    tmpdat$ha1gma0 <- ((I(tmpdat$m==1)*gma0 + I(tmpdat$m==0)*(1-gma0))/psm) * psa1 
    
    # targeting step using weights and logit(q1) as an offset
    # epsilon is calculated using full data, then applied to estimate of y under m=0 and m=1
    epsilona1gma0 <- coef(glm(y ~  1 , weights=ha1gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
    tmpdat$qyupm0a1gma0 <- plogis(qlogis(tmpdat$qyinit[,2]) + epsilona1gma0)
    tmpdat$qyupm1a1gma0 <- plogis(qlogis(tmpdat$qyinit[,3]) + epsilona1gma0)
    
    # marginalize over m 
    tmpdat$Qa1gma0 <- tmpdat$qyupm0a1gma0*(1-gma0) + tmpdat$qyupm1a1gma0*gma0
    
    # regress Q on covariates among those with a=1
    Qa1g0fit <- glm(formula = paste("Qa1gma0",qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
    Qa1g0 <- predict(Qa1g0fit, type="response", newdata=tmpdat)
    
    # update again (necessary when a is nonrandom)
    epsilonza1gma0 <-coef(glm(Qa1gma0~ 1 , weights=psa1, offset=qlogis(Qa1g0), family="quasibinomial", data=tmpdat))
    Qzupa1gma0 <-plogis(qlogis(Qa1g0) + epsilonza1gma0)
  
    # calculate tmle for y under a=1 and gm where a=0
    tmlea1m0 <- mean(Qzupa1gma0)

  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{0 , g_0})
  #-------------------------------------------------------------------------------------------------------------
    # calculate weights when a=0 and gm when a=0 
    tmpdat$ha0gma0 <- ((tmpdat$m*gma0 + (1-tmpdat$m)*(1-gma0))/psm) * psa0 
    
    # targeting step 
    epsilona0gma0 <- coef(glm(y ~  1 , weights=ha0gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
    tmpdat$qyupm0a0gma0 <- plogis(qlogis(tmpdat$qyinit[,2]) + epsilona0gma0)
    tmpdat$qyupm1a0gma0 <- plogis(qlogis(tmpdat$qyinit[,3]) + epsilona0gma0)
    
    # marginalize over m 
    tmpdat$Qa0gma0 <- tmpdat$qyupm0a0gma0*(1-gma0) + tmpdat$qyupm1a0gma0*gma0
    
    # regress Q on covariates among those with a=0
    Qa0g0fit <- glm(formula = paste("Qa0gma0",qmodel, sep="~"), data=tmpdat[tmpdat$a==0,], family="quasibinomial")
    Qa0g0 <- predict(Qa0g0fit, type="response", newdata=tmpdat)
    
    # update again (necessary when a is nonrandom)
    epsilonza0gma0 <- coef(glm(Qa0gma0~ 1 , weights=psa0, offset=qlogis(Qa0g0), family="quasibinomial", data=tmpdat))
    Qzupa0gma0 <-plogis(qlogis(Qa0g0) + epsilonza0gma0)
    
    
    # calculate tmle for y under a=0 and gm where a=0
    tmlea0m0 <- mean(Qzupa0gma0)

  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_1})
  #-------------------------------------------------------------------------------------------------------------
    # calculate weights when a=1 and gm when a=1 
    tmpdat$ha1gma1<- ((tmpdat$m*gma1 + (1-tmpdat$m)*(1-gma1))/psm) * psa1 
    
    # targeting step
    epsilona1gma1<-coef(glm(y ~  1 , weights=ha1gma1, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
    tmpdat$qyupm0a1gma1<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilona1gma1)
    tmpdat$qyupm1a1gma1<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilona1gma1)
    
    # marginalize over m
    tmpdat$Qa1gma1<-tmpdat$qyupm0a1gma1*(1-gma1) + tmpdat$qyupm1a1gma1*gma1
    
    # regress Q on covariates among those with a=1
    Qa1g1fit <- glm(formula = paste("Qa1gma1",qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
    Qa1g1 <- predict(Qa1g1fit, type="response", newdata=tmpdat)
    
    # update again (necessary when a is nonrandom)
    epsilonza1gma1 <-coef(glm(Qa1gma1~ 1 , weights=psa1, offset=qlogis(Qa1g1), family="quasibinomial", data=tmpdat))
    Qzupa1gma1 <- plogis(qlogis(Qa1g1) + epsilonza1gma1)
    
    # calculate tmle for y under a=1 and gm where a=1
    tmlea1m1 <- mean(Qzupa1gma1)

  #-------------------------------------------------------------------------------------------------------------
  #   estimate variances using EIC 
  #-------------------------------------------------------------------------------------------------------------
  
    tmpdat$qyupa1g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona1gma0)
    tmpdat$qyupa1g1<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona1gma1)
    tmpdat$qyupa0g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona0gma0)
    
    # EIC for E(Y_{1, g_0})
    eic1a1g0<-tmpdat$ha1gma0 * (tmpdat$y - tmpdat$qyupa1g0)
    eic2a1g0<-psa1*(tmpdat$Qa1gma0 - Qzupa1gma0)
    eic3a1g0<-Qzupa1gma0 -tmlea1m0
  
    eica1g0<-eic1a1g0 + eic2a1g0 + eic3a1g0
  
    # EIC for E(Y_{1, g_1})
    eic1a1g1<-tmpdat$ha1gma1 * (tmpdat$y - tmpdat$qyupa1g1)
    eic2a1g1<-psa1*(tmpdat$Qa1gma1- Qzupa1gma1)
    eic3a1g1<-Qzupa1gma1-tmlea1m1
  
    eica1g1<-eic1a1g1 + eic2a1g1 + eic3a1g1
  
    # EIC for E(Y_{0, g_0})
    eic1a0g0<-tmpdat$ha0gma0 * (tmpdat$y - tmpdat$qyupa0g0)
    eic2a0g0<-psa0*(tmpdat$Qa0gma0 - Qzupa0gma0)
    eic3a0g0<-Qzupa0gma0 -tmlea0m0
  
    eica0g0<-eic1a0g0 + eic2a0g0 + eic3a0g0
  
  
  #-------------------------------------------------------------------------------------------------------------
  #   combine to get mediation parameters, variances, and 95% confidence intervals  
  #-------------------------------------------------------------------------------------------------------------
  
    sde_tmle <- tmlea1m0-tmlea0m0
    sde_eic <- eica1g0 - eica0g0
    var_sde_eic <- var(sde_eic)/nrow(tmpdat)
    
    sie_tmle <- tmlea1m1-tmlea1m0
    sie_eic <- eica1g1 - eica1g0
    var_sie_eic <- var(sie_eic)/nrow(tmpdat)
    
    results <- data.frame(cbind(sde=sde_tmle, sde_var=var_sde_eic, sie=sie_tmle, sie_var=var_sie_eic, 
                                 sde_lb = sde_tmle - 1.96*sqrt(var_sde_eic), sde_ub = sde_tmle + 1.96*sqrt(var_sde_eic), 
                                 sie_lb = sie_tmle - 1.96*sqrt(var_sie_eic), sie_ub = sie_tmle + 1.96*sqrt(var_sie_eic)))
  
  return(results)
}
