library(pracma)
library(expm)
library(mvtnorm)
clear()
tic()
options(digits = 10)
w=2 #span MA(w=2)
n=5 #Sample size
u0=0 #In-control average
u1=0 #out-of-control average
u2=(u0+u1)/2
DESVIO <- 1/(n^0.5)
#The covariance matrix used in pmvnorm
Sigma2 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2), 
                   (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 2) 
Sigma3 <- matrix(c((DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2),   
                   0, (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w,(w-1)*(DESVIO^2)/(w^2),
                   0,(w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w), nrow = 3)
#We chose to leave the various conditional probabilities in their explicit form, 
#without using functions, for better understanding.
LSC_optim<- function(LSC){
  A <- Miwa()
  #A <- GenzBretz(abseps = tole, releps = 0)
  ####Branch 1#######
  mu2 <- c(u0, u0) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u0, u0, u2) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pa=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-LSC,LSC), upper = c(LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pas=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-LSC,-Inf), upper = c(LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pai=R3i[1] /R2[1]
  
  ####Branch 2#######
  mu2 <- c(u0, u2) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,LSC), upper = c(LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u0, u2, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,LSC,-LSC), upper=c(LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pb=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,LSC,LSC), upper=c(LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pbs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,LSC,-Inf), upper=c(LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pbi=R3i[1] /R2[1]
  
  
  ####Branch 3#######
  mu2 <- c(u0, u2) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-Inf), upper = c(LSC,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u0, u2, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-Inf,-LSC), upper=c(LSC,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pc=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-Inf,LSC), upper=c(LSC,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pcs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-Inf,-Inf), upper=c(LSC,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pci=R3i[1] /R2[1]
  
  ####Branch 4#######
  mu2 <- c(u0, u2) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u0, u2, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pd=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-LSC,LSC), upper=c(LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pds=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-LSC,-Inf), upper=c(LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pdi=R3i[1] /R2[1]
  
  ####Branch 5#######
  mu2 <- c(u2, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-Inf), upper = c(Inf,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-Inf,-LSC), upper=c(Inf,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pe=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-Inf,LSC), upper=c(Inf,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pes=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-Inf,-Inf), upper=c(Inf,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pei=R3i[1] /R2[1]
  
  ####Branch 7#######
  mu2 <- c(u2, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-LSC), upper = c(Inf,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-LSC,-LSC), upper=c(Inf,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pf=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-LSC,LSC), upper=c(Inf,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pfs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-LSC,-Inf), upper=c(Inf,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pfi=R3i[1] /R2[1]
  
  ####Branch 7#######
  mu2 <- c(u2, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,LSC), upper = c(-LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,LSC,-LSC), upper=c(-LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pg=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,LSC,LSC), upper=c(-LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pgs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,LSC,-Inf), upper=c(-LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pgi=R3i[1] /R2[1]
  
  ####Branch 10#######
  mu2 <- c(u2, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,-LSC), upper = c(-LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,-LSC,-LSC), upper=c(-LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  ph=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,-LSC,LSC), upper=c(-LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  phs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,-LSC,-Inf), upper=c(-LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  phi=R3i[1] /R2[1]
  
  ####Branch 11#######
  mu2 <- c(u2, u1)  # mean vector
  R2 <- pmvnorm(lower = c(-LSC,LSC), upper = c(LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # vetor média
  R3<-pmvnorm(lower = c(-LSC,LSC,-LSC), upper=c(LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pi=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,LSC,LSC), upper=c(LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pis=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,LSC,-Inf), upper=c(LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pii=R3i[1] /R2[1]
  
  ####Branch 12#######
  mu2 <- c(u2, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-Inf), upper = c(LSC,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-Inf,-LSC), upper=c(LSC,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pj=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-Inf,LSC), upper=c(LSC,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pjs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-Inf,-Inf), upper=c(LSC,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pji=R3i[1] /R2[1]
  
  
  ####Branch 13#######
  mu2 <- c(u2, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u2, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pk=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-LSC,LSC), upper=c(LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pks=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-LSC,-Inf), upper=c(LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pki=R3i[1] /R2[1]
  
  
  ####Branch 14#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,LSC), upper = c(LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,LSC,-LSC), upper=c(LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pl=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,LSC,LSC), upper=c(LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pls=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,LSC,-Inf), upper=c(LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pli=R3i[1] /R2[1]
  
  ####Branch 14#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-Inf), upper = c(LSC,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-Inf,-LSC), upper=c(LSC,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pm=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-Inf,LSC), upper=c(LSC,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pms=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-Inf,-Inf), upper=c(LSC,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pmi=R3i[1] /R2[1]
  
  
  ####Branch 16#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pn=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-LSC,LSC), upper=c(LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pns=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-LSC,-Inf), upper=c(LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pni=R3i[1] /R2[1]
  
  ####Branch 18#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-Inf), upper = c(Inf,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-Inf,-LSC), upper=c(Inf,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  po=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-Inf,LSC), upper=c(Inf,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pos=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-Inf,-Inf), upper=c(Inf,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  poi=R3i[1] /R2[1]
  
  ####Branch 19#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-LSC), upper = c(Inf,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-LSC,-LSC), upper=c(Inf,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pp=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-LSC,LSC), upper=c(Inf,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pps=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-LSC,-Inf), upper=c(Inf,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  ppi=R3i[1] /R2[1]
  
  ####Branch 20#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,LSC), upper = c(-LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1)# mean vector
  R3<-pmvnorm(lower = c(-Inf,LSC,-LSC), upper=c(-LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pq=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,LSC,LSC), upper=c(-LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pqs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,LSC,-Inf), upper=c(-LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pqi=R3i[1] /R2[1]
  
  
  ####Branch 22#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,-LSC), upper = c(-LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,-LSC,-LSC), upper=c(-LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  ps=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,-LSC,LSC), upper=c(-LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pss=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,-LSC,-Inf), upper=c(-LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  psi=R3i[1] /R2[1]
  
  
  ####Branch 23#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-Inf), upper = c(Inf,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-Inf,-LSC), upper=c(Inf,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pt=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-Inf,LSC), upper=c(Inf,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pts=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-Inf,-Inf), upper=c(Inf,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pti=R3i[1] /R2[1]
  
  
  ####Branch 25#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-LSC), upper = c(Inf,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-LSC,-LSC), upper=c(Inf,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pu=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-LSC,LSC), upper=c(Inf,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pus=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-LSC,-Inf), upper=c(Inf,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pui=R3i[1] /R2[1]
  
  
  ####Branch 26#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,LSC), upper = c(-LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,LSC,-LSC), upper=c(-LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pv=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,LSC,LSC), upper=c(-LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pvs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,LSC,-Inf), upper=c(-LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pvi=R3i[1] /R2[1]
  
  
  ####Branch 28#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,-LSC), upper = c(-LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,-LSC,-LSC), upper=c(-LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  px=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,-LSC,LSC), upper=c(-LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pxs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,-LSC,-Inf), upper=c(-LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pxi=R3i[1] /R2[1]
  
  
  ####Branch 29#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,LSC), upper = c(LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,LSC,-LSC), upper=c(LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pz=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,LSC,LSC), upper=c(LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pzs=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,LSC,-Inf), upper=c(LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pzi=R3i[1] /R2[1]
  
  ####Branch 30#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-Inf), upper = c(LSC,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-Inf,-LSC), upper=c(LSC,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pw=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-Inf,LSC), upper=c(LSC,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pws=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-Inf,-Inf), upper=c(LSC,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pwi=R3i[1] /R2[1]
  
  ####Branch 31#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-LSC,LSC), upper=c(LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  pys=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-LSC,-Inf), upper=c(LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  pyi=R3i[1] /R2[1]
  
  ####Branch 33#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-Inf), upper = c(Inf,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-Inf,-LSC), upper=c(Inf,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py1=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-Inf,LSC), upper=c(Inf,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py1s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-Inf,-Inf), upper=c(Inf,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py1i=R3i[1] /R2[1]
  
  
  ####Branch 34#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(LSC,-LSC), upper = c(Inf,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(LSC,-LSC,-LSC), upper=c(Inf,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py2=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(LSC,-LSC,LSC), upper=c(Inf,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py2s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(LSC,-LSC,-Inf), upper=c(Inf,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py2i=R3i[1] /R2[1]
  
  ####Branch 34#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,LSC), upper = c(-LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,LSC,-LSC), upper=c(-LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py3=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,LSC,LSC), upper=c(-LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py3s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,LSC,-Inf), upper=c(-LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py3i=R3i[1] /R2[1]
  
  
  ####Branch 37#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-Inf,-LSC), upper = c(-LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-Inf,-LSC,-LSC), upper=c(-LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py4=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-Inf,-LSC,LSC), upper=c(-LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py4s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-Inf,-LSC,-Inf), upper=c(-LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py4i=R3i[1] /R2[1]
  
  
  ####Branch 38#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,LSC), upper = c(LSC,Inf), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,LSC,-LSC), upper=c(LSC,Inf,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py5=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,LSC,LSC), upper=c(LSC,Inf,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py5s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,LSC,-Inf), upper=c(LSC,Inf,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py5i=R3i[1] /R2[1]
  
  
  ####Branch 39#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-Inf), upper = c(LSC,-LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-Inf,-LSC), upper=c(LSC,-LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py6=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-Inf,LSC), upper=c(LSC,-LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py6s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-Inf,-Inf), upper=c(LSC,-LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py6i=R3i[1] /R2[1]
  
  
  ####Branch 40#######
  mu2 <- c(u1, u1) # mean vector
  R2 <- pmvnorm(lower = c(-LSC,-LSC), upper = c(LSC,LSC), mean = mu2, sigma = Sigma2,algorithm = A)
  mu3 <- c(u1, u1, u1) # mean vector
  R3<-pmvnorm(lower = c(-LSC,-LSC,-LSC), upper=c(LSC,LSC,LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py7=R3[1] /R2[1]
  R3s<-pmvnorm(lower = c(-LSC,-LSC,LSC), upper=c(LSC,LSC,Inf), mean = mu3, sigma = Sigma3,algorithm = A)
  py7s=R3s[1] /R2[1]
  R3i<-pmvnorm(lower = c(-LSC,-LSC,-Inf), upper=c(LSC,LSC,-LSC), mean = mu3, sigma = Sigma3,algorithm = A)
  py7i=R3i[1] /R2[1]
  
  # Chain Markov
  c1 <- c(0, pas, pai, pa, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c2 <- c(0, 0, 0, 0, pbs, pbi, pb, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c3 <- c(0, 0, 0, 0, 0, 0, 0, pcs, pci, pc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c4 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pds, pdi, pd, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c5 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c6 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pes, pei, pe, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c7 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pfs, pfi, pf, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c8 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pgs, pgi, pg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c9 <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c10 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, phs, phi, ph, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c11 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pis, pii, pi, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c12 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pjs, pji, pj, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c13 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pks, pki, pk, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c14 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pls, pli, pl, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c15 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pms, pmi, pm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c16 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pns, pni, pn, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c17 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c18 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pos, poi, po, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c19 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pps, ppi, pp, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c20 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pqs, pqi, pq, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c21 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c22 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pss, psi, ps, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c23 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c24 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pts, pti, pt, 0, 0, 0, 0)
  c25 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, pus, pui, pu)
  c26 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pvs, pvi, pv, 0, 0, 0, 0, 0, 0)
  c27 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  c28 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pxs, pxi, px)
  c29 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pzs, pzi, pz, 0, 0, 0, 0, 0, 0)
  c30 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pws, pwi, pw, 0, 0, 0)
  c31 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pys, pyi, py)
  c32 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
  c33 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py1s, py1i, py1, 0, 0, 0)
  c34 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py2s, py2i, py2)
  c35 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py3s, py3i, py3, 0, 0, 0, 0, 0, 0)
  c36 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
  c37 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py4s, py4i, py4)
  c38 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py5s, py5i, py5, 0, 0, 0, 0, 0, 0)
  c39 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py6s, py6i, py6, 0, 0, 0)
  c40 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, py7s, py7i, py7)
  
  
  M <- rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
             c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c32, c33, c34, c35, c36, c37, c38, c39, c40)
  M <- M / rowSums(M)
  
  #The transient submatrix Q of the Markov chain
  A <- c(5, 9, 17, 21, 23, 27,32,36) #Absorbing states
  Q <- M[-A, -A]
  
  #Obtaining the ARL
  V=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  I=diag(rep(1,32)) 
  U=ones(32,1)
  P1=solve((I-Q))
  ARLa=V%*%P1
  ARLb=ARLa%*%U
  return((ARLb-370.4)^2)
  
}

LI=0.5
LS=2.5

par_optim<- optimize(LSC_optim,lower=LI,upper = LS,tol = .Machine$double.eps)

Ua=par_optim[[1]][1] 
cat('UCL=',Ua,"\n") 
cat('LCL=',-Ua,"\n") 
toc()
