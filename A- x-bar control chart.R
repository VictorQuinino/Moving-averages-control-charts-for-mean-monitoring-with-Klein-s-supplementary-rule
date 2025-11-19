#A control chart for the monitoring of a process mean
rm(list = ls())

u0=0#In-control average
u1=0 #out-of-control average 
s0=1 #in-control standard deviation 
n=5  #Sample size
ARL0=370.4 #Target ARL0 

alfaa=1/ARL0 #Prob Xbar

LSCxb=qnorm((1-alfaa/2),u0,s0/(n^0.5)) #upper control limit
LICxb=qnorm(alfaa/2,u0,s0/(n^0.5)) # lower control limit

#out-of-control
pa=pnorm(LICxb,u1,s0/(n^0.5)) #probability of being below LCL
pb=1-pnorm(LSCxb,u1,s0/(n^0.5)) #probability of being above UCL
pc=1-pa-pb #probability of being between control limits
pd=pa+pb
ARL1<- 1/(pd)

#SDRL
SDRL=((1-pd)/(pd^2))^0.5


cat('LSCxb=',LSCxb,"\n")
cat('LICxb=',LICxb,"\n")
cat('Probxbar=',alfaa,"\n")
cat('ARL1=',ARL1,"\n")
cat('SDRL=',SDRL,"\n")