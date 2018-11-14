mle.logreg.constrained = function(fmla, data) {
    # Define the negative log likelihood function
    logl <- function(theta,x,y){
      y <- y
      x <- as.matrix(x)
      beta <- theta[1:ncol(x)]

      # Use the log-likelihood of the Bernouilli distribution, where p is
      # defined as the logistic transformation of a linear combination
      # of predictors, according to logit(p)=(x%*%beta)
      loglik <- sum(-y*log(1 + exp(-(x%*%beta)))-(1-y)*log(1 + exp(x%*%beta)))
      return(-loglik)
    }

    # Prepare the data
    outcome = rownames(attr(terms(fmla),"factors"))[1]
    dfrTmp = model.frame(data)
    x = as.matrix(model.matrix(fmla, data=dfrTmp))
    y = as.matrix(data[,match(outcome,colnames(data))])

    # Define initial values for the parameters
    theta.start = rep(0,(dim(x)[2]))
    names(theta.start) = colnames(x)
    # print("theta.start:"); print(theta.start)

    # Non-negative slopes constraint
    ## ************************************************************************
    ## OS: MODIFIED TO KARA'S PROBLEM TO RESTRICT \beta_1 >= 0.0001
    ## ALL OTHER COEFFICIENTS (including intercept) are UNCONSTRAINED
    lower = c(-Inf,0.0001,rep(-Inf,(length(theta.start)-2)))
    ## ************************************************************************

    # Calculate the maximum likelihood
    mle = optim(theta.start,logl,x=x,y=y,hessian=T,lower=lower,
      method="L-BFGS-B")

    # Obtain regression coefficients
    beta = mle$par

    # Calculate the Information matrix
    # The variance of a Bernouilli distribution is given by p(1-p)
    p = 1/(1+exp(-x%*%beta))
    V = array(0,dim=c(dim(x)[1],dim(x)[1]))
    diag(V) = p*(1-p)
    IB = t(x)%*%V%*%x

    # Return estimates
    out = list(beta=beta,vcov=solve(IB),dev=2*mle$value)
}

predictClogis = function(newX, beta) {
 eta <- cbind(1,as.matrix(newX))[,!is.na(beta),
 drop = FALSE] %*% beta[!is.na(beta)]
 as.vector(plogis(eta))
}
