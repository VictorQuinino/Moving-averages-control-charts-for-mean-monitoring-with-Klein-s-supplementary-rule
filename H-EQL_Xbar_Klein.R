library(cubature)
library(pracma)
tic()
clear ()
s = 1 #in-control standard deviation
n = 1 #Sample size
#Klein Limits
LSC=1.781418898 #n=1
#LSC=0.7966747505 #n=5
LIC=-LSC
#Limits for Double Integrals
lowerLimit <- c(0)  # Lower bound for u 
upperLimit <- c(0.75)  # Upper bound for u 
#Uniform Distribution
R1=(1/(upperLimit[1]-lowerLimit[1]))
# Function that calculates the EQL 
Inte <- function(params) {
u1=params[1]
pxi=pnorm(LIC,u1,s/(n^0.5))
pxs=1-pnorm(LSC,u1,s/(n^0.5))
pxc=1-pxi-pxs
  
size<- 5 #Size of the markov chain
MarkovChain<- matrix(0,nrow=size,ncol=size,byrow=TRUE)
MarkovChain[1,]<- c(pxc,pxs,pxi,0,0) 
MarkovChain[2,]<- c(pxc,0,pxi,pxs,0) 
MarkovChain[3,]<- c(pxc,pxs,0,0,pxi) 
MarkovChain[4,]<- c(0,0,0,1,0) 
MarkovChain[5,]<- c(0,0,0,0,1) 
#solving the markov chain through absorbing states
n_states <- nrow(MarkovChain) 
transient_idx <- 1:3 # xc, xs, xi 
absorbing_idx <- 4:5 # xss, xii
R <- MarkovChain[transient_idx, transient_idx, drop = FALSE]
u <- c(1,0,0)
I<-diag(rep(1,3)) 
ones<-ones(3,1)
P1<-solve((I-R))
ARL1<-u%*%P1%*%ones
P2=(u1^2)*ARL1
return(P2)
}
EQL_Result <- adaptIntegrate(Inte, lowerLimit = lowerLimit, upperLimit = upperLimit)
EQLA=EQL_Result$integral*(R1)
cat("EQL_Klein",EQLA,"\n")
toc()