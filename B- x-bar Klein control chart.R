#A control chart for the monitoring of mean with Klein's supplementary rule
library(pracma)
library(expm)
rm(list = ls())
####################################################################################################
u0=0 #In-control average
u1=0 #out-of-control average 
s0=1 #in-control standard deviation 
s1=1 # out-of-control standard deviation 
n=5  #Sample size
ARL0= 370.4 #Target ARL0

####################################################################################################
OtiUCL <- function(U){#Optimization function used to find the UCL and LCL
  UCLxb=qnorm((1-U/2),u0,s0/(n^0.5))
  LCLxb=qnorm(U/2,u0,s0/(n^0.5))
  #X-bar control chart 
  pxi=pnorm(LCLxb,u0,s0/(n^0.5))
  pxs=1-pnorm(UCLxb,u0,s0/(n^0.5))
  pxc=1-pxs-pxi
  
  #Markov chain
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
  I=diag(rep(1,3)) 
  ones=ones(3,1)
  P1=solve((I-R))
  
  ARL=u%*%P1%*%ones
  
  ARLphi=(ARL-ARL0)^2
  return(ARLphi)
}

####################################################################################################
#Limit used in the optimize function for the x-bar control chart 
LSa=1 
#initial value. Change as needed.
ivxbar=0.04
par_optim <- nlminb(ivxbar,OtiUCL,lower=1e-6,upper =LSa)#

Ua <- par_optim$par

UCLxb=qnorm((1-Ua/2),u0,s0/(n^0.5)) #upper control limit for the x-bar
LCLxb=qnorm(Ua/2,u0,s0/(n^0.5)) #lower control limit for the x-bar

#X-bar control chart
pxi=pnorm(LCLxb,u1,s1/(n^0.5))
pxs=1-pnorm(UCLxb,u1,s1/(n^0.5))
pxc=1-pxs-pxi

####################################################################################################
#Now that we have the probabilities for this specific case, we solve it through a Markov Chain again
#Markov chain
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

ARL<-u%*%P1%*%ones

M1<-P1%^%2

SDRL<-(2%*%u%*%M1%*%R%*%ones-(ARL^2)+ARL)^0.5

###############################################################

options(digits=7)
cat('UCLxb=',UCLxb,"\n")
cat('LCLxb=',LCLxb,"\n")
cat('ARL=',ARL,"\n")
cat('SDRL=',SDRL,"\n")