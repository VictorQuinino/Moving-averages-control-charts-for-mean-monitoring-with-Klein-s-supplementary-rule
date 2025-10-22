#A control chart for the monitoring of a process mean
rm(list = ls())

u0=0#In-control average
u1=0 #out-of-control average 
s0=1 #in-control standard deviation 
s1=1 # out-of-control standard deviation
n=5  #Sample size
ARL0=370.4 #Target ARL0 


OtiUCL <- function(alfa){#simple optimization function so that alpha always leads to desired ARL0
  ARL=1/alfa
  ARLphi=(ARL-ARL0)^2
  return(ARLphi)
}
#upper limit to be used in the optimize function. Change as needed.
LSa=1 
LSb=1
#initial value for optimization, Change as needed.
ivxbar=0.0442
ivs2=0.0442
par_optim <- nlminb(c(ivxbar,ivs2),OtiUCL,lower=c(0,0),upper =c(LSa,LSb))
alfaa=par_optim[[1]][1] #Prob Xbar

LSCxb=qnorm((1-alfaa/2),u0,s0/(n^0.5)) #upper control limit
LICxb=qnorm(alfaa/2,u0,s0/(n^0.5)) # lower control limit

pc=pnorm(LICxb,u1,s1/(n^0.5)) #probability of being below LCL
pa=1-pnorm(LSCxb,u1,s1/(n^0.5)) #probability of being above UCL
pb=1-pa-pc #probability of being between control limits


ARL1<- 1/(1-pb)

#SDRL
SDRL=((pb)^0.5)/(1-pb)


cat('LSCxb=',LSCxb,"\n")
cat('LICxb=',LICxb,"\n")
cat('Probxbar=',alfaa,"\n")
cat('ARL1=',ARL1,"\n")
cat('SDRL=',SDRL,"\n")