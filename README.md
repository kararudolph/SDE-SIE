# SDE-SIE
Estimators for stochastic direct and indirect effects.

Consider observed data O=(W,A,Z,M,Y) for each of n individuals, where we assume O<sub>1</sub>, ..., O<sub>n</sub>, are i.i.d. for the 
true, unknown data distribution, P<sub>0</sub>, on O. The subscript 0 denotes values under this true, unknown distribution P<sub>0</sub>. 
P is any probability distribution in statistical model, which is the set of distributions for which our estimand is 
identifiable. 

* W is a vector of exogenous baseline covariates, W=f(U<sub>W</sub>), where U<sub>W</sub> is unobserved exogenous error on W.
* A is a binary treatment, which may be an instrumental variable for several of the functions A=f(W, U<sub>A</sub>). 
* Z is an intermediate variable, Z=f(W, A, U<sub>Z</sub>).
* M is a mediator. Depending on the function, M=f(W, Z, U<sub>M</sub>) or M=f(W, A, Z, U<sub>M</sub>), see details in each of the functions.  
* Y is an outcome. Depending on the function, Y=f(W,Z,M,U<sub>Y</sub>) or Y=f(W,A,Z,M,U<sub>Y</sub>), see details in each of the functions. 

ivmedtmle.R estimates the 1) stochastic direct effect of A on Y not through M and the 2) indirect effect of A to Z to M to Y 
in the scenario where A is an instrumental variable which adheres to exclusion restriction M=f(W, Z, U<sub>M</sub>) and Y=f(W,Z,M,U<sub>Y</sub>).

medtmle_intermedvar.R estimates the 1) stochastic direct effect of A on Y not through M and the 2) indirect effect of A 
to Z to M to Y in the scenario where A is not an instrumental variable, so there is no exclusion restriction: M=f(W, A, Z, U<sub>M</sub>) and 
Y=f(W,A,Z,M,U<sub>Y</sub>).

medtmle_nointermedvar.R estimates the 1) stochastic direct effect of A on Y not through M and the 2) indirect effect of 
A to M to Y in the scenario where A is not an instrumental variable and there is no intermediate variable, Z,
so there is no exclusion restriction: M=f(W, A, U<sub>M</sub>) and Y=f(W,A,M,U<sub>Y</sub>).

CSDE_*.R functions estimate the complier stochastic direct effect of Z on Y not through M, using A as an instrument to 
address observed and unobserved confounding in the scenario where A is an A is an instrumental variable which adheres 
to exclusion restriction M=f(W, Z, U<sub>M</sub>) and Y=f(W,Z,M,U<sub>Y</sub>).
