library(pracma)
library(expm)
library(mvtnorm)
library(writexl)
clear()
tic()
tole=1e-8
options(digits = 5)
#A <- GenzBretz(abseps = tole, releps = 0)
A <- Miwa()
w=2 #span MA(w=2)
n=5 #Sample size
u0=0 #In-control average
#Use Program C
LSC=0.9428364267 #n=5
#LSC=2.108246345  #n=1

DESVIO <- 1/(n^0.5) #in-control standard deviation
#u1 out-of-control average 
UU=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5)
T=length(UU)
Result=matrix(0,T,3)
Sigma2 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2), 
                   (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 2) # matriz de covariância
Sigma3 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2),   
                   0, (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w,(w-1)*(DESVIO^2)/(w^2),
                   0,(w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 3)
for (i in 1:T){
  A <- Miwa()
  u1=UU[i]
  u2=(u0+u1)/2
  DESVIO <- 1/(n^0.5)
  
  ####Branch 1#######
  mu2 <- c(u0, u0) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, 
                sigma = Sigma2,algorithm = A)
  mu3 <- c(u0, u0, u2) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, 
              sigma = Sigma3,algorithm = A)
  p1c=R3[1] /R2[1]
  ####Branch 2#######
  mu2a <- c(u0, u2) # mean vector
  R2a <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2a, 
                 sigma = Sigma2,algorithm = A)
  mu3a <- c(u0, u2, u1) # mean vector
  R3a<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3a, 
               sigma = Sigma3,algorithm = A)
  p2c=R3a[1] /R2a[1]
  
  ####Branch 3#######
  mu2b <- c(u2, u1) # mean vector
  R2b <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2b, 
                 sigma = Sigma2,algorithm = A)
  mu3b <- c(u2, u1, u1) # mean vector
  R3b<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3b, 
               sigma = Sigma3,algorithm = A)
  p3c=R3b[1] /R2b[1]
  
  ####Branch 4#######
  mu2c <- c(u1, u1) # mean vector
  R2c <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2c, 
                 sigma = Sigma2,algorithm = A)
  mu3c <- c(u1, u1, u1) # mean vector
  R3c<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3c, 
               sigma = Sigma3,algorithm = A)
  p4c=R3c[1] /R2c[1]
  
  #Chain Markov
  c1<-c(0,p1c,(1-p1c),0,0,0,0,0,0)
  c2<-c(0,0,0,p2c,(1-p2c),0,0,0,0)
  c3<-c(0,0,1,0,0,0,0,0,0)
  c4<-c(0,0,0,0,0,p3c,(1-p3c),0,0)
  c5<-c(0,0,0,0,1,0,0,0,0)
  c6<-c(0,0,0,0,0,0,0,p4c,(1-p4c))
  c7<-c(0,0,0,0,0,0,1,0,0)
  c8<-c(0,0,0,0,0,0,0,p4c,(1-p4c))
  c9<-c(0,0,0,0,0,0,0,0,1)
  M<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9)
  
  #The transient submatrix Q of the Markov chain
  A=c(3,5,7,9) #Absorbing states
  Q<-M[-A,-A]
  #Obtaining the ARL
  V <- c(1,0,0,0,0)
  I=diag(rep(1,5)) 
  U=ones(5,1)
  P1=solve((I-Q))
  ARLa=V%*%P1
  ARLb=ARLa%*%U
  
  AR=V%*%P1%*%U
  M1=P1%^%2
  SDRL=(2%*%V%*%M1%*%Q%*%U-(AR^2)+AR)^0.5
  
  Result[i,1]=u1
  Result[i,2]=ARLb
  Result[i,3]=SDRL
  
}

Result=data.frame(Result)
writexl::write_xlsx(Result,"Result.xlsx")

toc()

