library(pracma)
library(mvtnorm)
clear()
tic()
tole=1e-8
options(digits = 10)
#A <- GenzBretz(abseps = tole, releps = 0)
Alvo=370.4
A <- Miwa()
w=2 #span MA(w=2)
n=5  #Sample size
u0=0 #In-control average
u1=0 #out-of-control average 
u2=(u1+u0)/2
DESVIO <- 1/(n^0.5)
#The covariance matrix used in pmvnorm
Sigma2 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2), 
                   (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 2) 
Sigma3 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2),   
                   0, (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w,(w-1)*(DESVIO^2)/(w^2),
                   0,(w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 3)
LSC_optim<- function(LSC){
  A <- Miwa()
  ####Branch 1#######
  mu2 <- c(u0, u0) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, 
                sigma = Sigma2,algorithm = A)
  mu3 <- c(u0, u0, u2) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, 
              sigma = Sigma3,algorithm = A)
  p1c=R3[1] /R2[1]
  
  #####Branch 2#######
  mu2a <- c(u0, u2) # mean vector
  R2a <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2a, 
                 sigma = Sigma2,algorithm = A)
  mu3a <- c(u0, u2, u1) # mean vector
  R3a<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3a, 
               sigma = Sigma3,algorithm = A)
  p2c=R3a[1] /R2a[1]
  
  #####Branch 3#######
  
  mu2b <- c(u2, u1) # mean vector
  R2b <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2b, 
                 sigma = Sigma2,algorithm = A)
  mu3b <- c(u2, u1, u1) # mean vector
  R3b<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3b, 
               sigma = Sigma3,algorithm = A)
  p3c=R3b[1] /R2b[1]
  
  #####Branch 4#######
  mu2c <- c(u1, u1) # mean vector
  R2c <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2c, 
                 sigma = Sigma2,algorithm = A)
  mu3c <- c(u1, u1, u1) # vetor média
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
  return((ARLb-Alvo)^2)
}
#Range for the control limits. Modify if necessary
LI=0.5
LS=2.5

par_optim<- optimize(LSC_optim,lower=LI,upper = LS,tol = .Machine$double.eps)

Ua=par_optim[[1]][1] 
cat('UCL=',Ua,"\n") 
cat('LCL=',-Ua,"\n") 
toc()

      