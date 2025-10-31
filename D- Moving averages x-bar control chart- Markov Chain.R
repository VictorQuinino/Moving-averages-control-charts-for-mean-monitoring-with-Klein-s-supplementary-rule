library(pracma)
library(expm)
library(mvtnorm)
library(cubature)
library(writexl)
clear()
tic()
tole=1e-8
options(digits = 5)
#A <- GenzBretz(abseps = tole, releps = 0)
A <- Miwa()
w=2
n=1
u0=0
#LSC=0.9428364276  #n=5
LSC=2.108246345  #n=1
DESVIO <- 1/(n^0.5)

UU=c(0,0.1,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5)

Result=matrix(0,10,3)
Sigma2 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2), 
(w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 2) # matriz de covariância
Sigma3 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2),   
                   0, (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w,(w-1)*(DESVIO^2)/(w^2),
                   0,(w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 3)
for (i in 1:10){
A <- Miwa()
u1=UU[i]
u2=(u0+u1)/2
DESVIO <- 1/(n^0.5)

mu2 <- c(u0, u0) # vetor média

R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)

mu3 <- c(u0, u0, u2) # vetor média
R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)

pa=R3[1] /R2[1]


###################Ramo 2##################

mu2a <- c(u0, u2) # vetor média

R2a <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2a, sigma = Sigma2,algorithm = A)


mu3a <- c(u0, u2, u1) # vetor média

R3a<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3a, sigma = Sigma3,algorithm = A)


pb=R3a[1] /R2a[1]

###################Ramo 3##################

mu2b <- c(u2, u1) # vetor média

R2b <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2b, sigma = Sigma2,algorithm = A)


mu3b <- c(u2, u1, u1) # vetor média
R3b<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3b, sigma = Sigma3,algorithm = A)


pc=R3b[1] /R2b[1]

####Ramo4#####

mu2c <- c(u1, u1) # vetor média

R2c <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2c, sigma = Sigma2,algorithm = A)

mu3c <- c(u1, u1, u1) # vetor média
R3c<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3c, sigma = Sigma3,algorithm = A)

pd=R3c[1] /R2c[1]

c7 <- c(0, pa, 0,0,0)   
c8 <- c(0, 0,pb,0,0)   
c9 <- c(0, 0,0,pc,0)
c10<-c(0, 0,0,0,pd)
c11<-c(0, 0,0,0,pd)
Q  <- rbind(c7, c8, c9,c10,c11)

V <- c(1,0,0,0,0)
I=diag(rep(1,5)) 
U=ones(5,1)
P1=solve((I-Q))
ARLa=V%*%P1
ARLb=ARLa%*%U


AR=V%*%P1%*%U
M1=P1%^%2
SDRL=(2%*%V%*%M1%*%Q%*%U-(AR^2)+AR)^0.5





#cat('ARL_Markov=',ARLb,"\n") 
Result[i,1]=u1
Result[i,2]=ARLb
Result[i,3]=SDRL
}

Result=data.frame(Result)
writexl::write_xlsx(Result,"Result.xlsx")

toc()

      