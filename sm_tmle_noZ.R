#-----------------------------------------------------------------------------------------------------------------------------
#     
# This is a function to estimate stochastic direct and indirect effects when there are direct effects of the 
#   exposure on the mediator and the outcome, and there is no mediator-outcome confounder affected by prior exposure.

# The stochastic direct effect is defined as E(Y_{1,g_0}) - E(Y_{0,g_0})
#   and the stochastic indirect effect is defined as E(Y_{1,g_1}) - E(Y_{1,g_0})
#   where Y_{a,g_a*} is the counterfactual value of Y when a=a and the mediator m is drawn from 
#   the distribution g under a*. 
#
# Inputs to the function: 
#     1. obsdat 
#       - should be a data frame that includes an exposure variable (called a with values 0/1), 
#              a mediator variable (called m with values 0/1), and an outcome variable (called y with values 0/1). 
#       - there can be other covariates that can be named anything and have any values.
#     2. amodel 
#       - the parametric model for a that includes the parents of the exposure. For example: "a ~ w1 + w2 + w3"
#     3. mmodel 
#       - the parametric model for m that includes the parents of the mediator. 
#     4. ymodel 
#       - the parametric model for y that includes the parents of the outcome. 
#     5. qmodel 
#       - the covariates to marginalize over. For example, "w1 + w2 + w3"
#  
#      
#       
#-----------------------------------------------------------------------------------------------------------------------------

sm_tmle_noZ <- function(obsdat, amodel, mmodel, ymodel, qmodel) {

  dfa1 <- dfa0 <- dfm1 <- dfm0 <- tmpdat <- obsdat
  
  dfa1$a <- dfm1$m <-  1 
  dfa0$a <- dfm0$m <-  0 

  
  # calculate gm under a=0 and gm under a=1
  gma1 <- predict(glm(formula=mmodel, family="binomial", data=tmpdat), newdata=dfa1, type="response")
  gma0 <- predict(glm(formula=mmodel, family="binomial", data=tmpdat), newdata=dfa0, type="response")
  
  # calculate a weights 
  psa1 <- I(tmpdat$a==1)/predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response")
  psa0 <- I(tmpdat$a==0)/(1-predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response"))
  
  # calculate m weights 
  ma <- predict(glm(formula=mmodel, family="binomial", data=tmpdat), type="response")
  psm <- (ma*tmpdat$m) + ((1-ma)*(1-tmpdat$m))


  # first generate initial estimate of y in full pop, and then under m=0 and m=1 
  tmpdat$qyinit <- cbind(predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=tmpdat, type="response"), 
                        predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm0, type="response"),
                        predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm1, type="response"))
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_0})
  #-------------------------------------------------------------------------------------------------------------
    # calculate weights when a=1 and gm when a=0 
    tmpdat$ha1gma0 <- ((tmpdat$m*gma0 + (1-tmpdat$m)*(1-gma0))/psm) * psa1 
  
    # targeting step using weights and logit(q1) as an offset
    # epsilon is calculated using full data, then applied to estimate of y under m=0 and m=1
    epsilona1gma0 <-coef(glm(y ~  1 , weights=ha1gma0 , offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
    tmpdat$qyupm0a1gma0 <- plogis(qlogis(tmpdat$qyinit[,2]) + epsilona1gma0)
    tmpdat$qyupm1a1gma0 <- plogis(qlogis(tmpdat$qyinit[,3]) + epsilona1gma0)
    
    # marginalize over m 
    tmpdat$Qa1gma0 <-tmpdat$qyupm0a1gma0*(1-gma0) + tmpdat$qyupm1a1gma0*gma0
    
    # regress Q on W among those with A=1 and predict
    Qzfita1gma0 <- glm(formula=paste("Qa1gma0", qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
    Qza1gma0 <-predict(Qzfita1gma0, type="response", newdata=tmpdat)
    
    # update again (necessary when a is nonrandom)
    epsilonza1gma0 <- coef(glm(Qa1gma0 ~ 1 , weights=psa1, offset=qlogis(Qza1gma0), family="quasibinomial", data=tmpdat))
    Qzupa1gma0 <- plogis(qlogis(Qza1gma0) + epsilonza1gma0)
    
    # calculate tmle for y under a=1 and gm where a=0
    tmlea1m0 <- mean(Qzupa1gma0)
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{0 , g_0})
  #-------------------------------------------------------------------------------------------------------------
    # calculate weights when a=0 and gm when a=0 
    tmpdat$ha0gma0 <- ((tmpdat$m*gma0 + (1-tmpdat$m)*(1-gma0))/psm) * psa0 
    
    # targeting step 
    epsilona0gma0 <-coef(glm(y ~  1 , weights=ha0gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
    tmpdat$qyupm0a0gma0 <-plogis(qlogis(tmpdat$qyinit[,2]) + epsilona0gma0)
    tmpdat$qyupm1a0gma0 <-plogis(qlogis(tmpdat$qyinit[,3]) + epsilona0gma0)
    
    # marginalize over m 
    tmpdat$Qa0gma0 <-tmpdat$qyupm0a0gma0*(1-gma0) + tmpdat$qyupm1a0gma0*gma0
    
    # regress Q on W among those with A=0 and predict
    Qzfita0gma0 <- glm(formula=paste("Qa0gma0", qmodel, sep="~"), data=tmpdat[tmpdat$a==0,], family="quasibinomial")
    Qza0gma0 <- predict(Qzfita0gma0, type="response", newdata=tmpdat)
    
    # update again (necessary when a is nonrandom)
    epsilonza0gma0 <- coef(glm(Qa0gma0 ~ 1 , weights=psa0, offset=qlogis(Qza0gma0), family="quasibinomial", data=tmpdat))
    Qzupa0gma0 <- plogis(qlogis(Qza0gma0) + epsilonza0gma0)
    
    # calculate tmle for y under a=0 and gm where a=0
    tmlea0m0 <- mean(Qzupa0gma0)
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_1})
  #-------------------------------------------------------------------------------------------------------------
    # calculate weights when a=1 and gm when a=1 
    tmpdat$ha1gma1 <- ((tmpdat$m*gma1 + (1-tmpdat$m)*(1-gma1))/psm) * psa1 
     
    # targeting step
    epsilonma1gma1<-coef(glm(y ~  1 , weights=ha1gma1, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
    tmpdat$qyupm0a1gma1<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilonma1gma1)
    tmpdat$qyupm1a1gma1<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilonma1gma1)
     
    # marginalize over m
    tmpdat$Qa1gma1<-tmpdat$qyupm0a1gma1*(1-gma1) + tmpdat$qyupm1a1gma1*gma1
     
    # regress Q on W among those with A=1 and predict
    Qzfita1gma1<-glm(formula=paste("Qa1gma1", qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
    Qza1gma1<-predict(Qzfita1gma1, type="response", newdata=tmpdat)
    
    # update again (necessary when a is nonrandom)
    epsilonza1gma1<-coef(glm(Qa1gma1~ 1 , weights=psa1, offset=qlogis(Qza1gma1), family="quasibinomial", data=tmpdat))
    Qzupa1gma1<-plogis(qlogis(Qza1gma1) + epsilonza1gma1)
    
    # calculate tmle for y under a=1 and gm where a=1
    tmlea1m1 <- mean(Qzupa1gma1)
 
  #-------------------------------------------------------------------------------------------------------------
  #   estimate variances using EIC  
  #-------------------------------------------------------------------------------------------------------------

     tmpdat$qyupa1g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona1gma0)
     tmpdat$qyupa1g1<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilonma1gma1)
     tmpdat$qyupa0g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona0gma0)
     
     # EIC for E(Y_{1, g_0})
     eic1a1g0 <- tmpdat$ha1gma0 * (tmpdat$y - tmpdat$qyupa1g0)
     eic2a1g0 <- psa1*(tmpdat$Qa1gma0 - Qzupa1gma0)
     eic3a1g0 <- Qzupa1gma0 -tmlea1m0
     
     eica1g0 <- eic1a1g0 + eic2a1g0 + eic3a1g0
     
     # EIC for E(Y_{1, g_1})
     eic1a1g1 <- tmpdat$ha1gma1 * (tmpdat$y - tmpdat$qyupa1g1)
     eic2a1g1 <- psa1*(tmpdat$Qa1gma1- Qzupa1gma1)
     eic3a1g1 <- Qzupa1gma1-tmlea1m1
     
     eica1g1 <- eic1a1g1 + eic2a1g1 + eic3a1g1
     
     # EIC for E(Y_{0, g_0})
     eic1a0g0 <- tmpdat$ha0gma0 * (tmpdat$y - tmpdat$qyupa0g0)
     eic2a0g0 <- psa0*(tmpdat$Qa0gma0 - Qzupa0gma0)
     eic3a0g0 <- Qzupa0gma0 - tmlea0m0
     
     eica0g0 <- eic1a0g0 + eic2a0g0 + eic3a0g0
     
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
