library(cubature)
library(pracma)
tic()
clear ()
s0 = 1 #In-control
u0 = 0 #In-control
n = 1 #Sample size
ARL0=370.4
#Limits
alfaa=1/ARL0 #Prob Xbar
LSC=qnorm((1-alfaa/2),u0,s0/(n^0.5)) #upper control limit
LIC=qnorm(alfaa/2,u0,s0/(n^0.5)) # lower control limit
#Limits for Integral
lowerLimit <- c(0)  # Lower bound for u 
upperLimit <- c(0.75)  # Upper bound for u 
#Uniform Distribution
R1=(1/(upperLimit[1]-lowerLimit[1]))
# Function that calculates the EQL 
Inte <- function(params) {
u1=params[1]
pi=pnorm(LIC,u1,s0/(n^0.5))
ps=1-pnorm(LSC,u1,s0/(n^0.5))
pc=1-pi-ps
ARL1=1/(pi+ps)
P2=(u1^2)*ARL1
return(P2)
}
EQL_Result <- adaptIntegrate(Inte, lowerLimit = lowerLimit, upperLimit = upperLimit)
EQLA=EQL_Result$integral*(R1)
cat("EQL_Klein",EQLA,"\n")
toc()