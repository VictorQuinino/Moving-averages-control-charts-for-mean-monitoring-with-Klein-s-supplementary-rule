library(mvtnorm)
library(pracma)
clear()
tic()
n_sim <- 100000   #number of Monte Carlo replications
u0 <- #In-control average
  n <- 1
w <- 2
DESVIO <- 1 / sqrt(n)
LSC <- 2.108246345 
LIC=-LSC
Sigma3 <- matrix(c(
  (DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2), 0,
  (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w, (w-1)*(DESVIO^2)/(w^2),
  0, (w-1)*(DESVIO^2)/(w^2), (DESVIO^2)/w
), nrow = 3, byrow = TRUE)
#out-of-control average
UU=c(0,0.2,0.4,0.6,0.8,1,1.25,1.5)
T=length(UU)
Result=matrix(0,T,3)
for (ii in 1:T){
  print(ii)
  u1=UU[ii]
  u2=(u0+u1)/2
  mu3  <- c(u0, u0, u2)
  mu3a <- c(u0, u2, u1)
  mu3b <- c(u2, u1, u1)
  mu3c <- c(u1, u1, u1)
  GG=matrix(0,n_sim,1)
  for (i in 1:n_sim){
    D<-c()
    s1<-0
    z=0
    while (z==0){
      sim_data <- rmvnorm(1, mean = mu3, sigma = Sigma3)
      
      controle_prev <- ((sim_data[,1] <= LSC & sim_data[,1] >= -LSC) &
                          (sim_data[,2] <= LSC & sim_data[,2] >= -LSC))
      z=z+controle_prev
    }
    
    controle_prevv <-(sim_data[,3] <= LSC & sim_data[,3] >= -LSC)
    
    if (controle_prevv==0){
      D=c(D,1)
      s1=s1+1
      GG[i]=sum(D)
      next
    }
    
    D=c(D,1)
    z1=0
    while (z1==0){
      sim_data1 <- rmvnorm(1, mean = mu3a, sigma = Sigma3)
      controle_prev1 <- ((sim_data1[,1] <= LSC & sim_data1[,1] >= -LSC) &
                           (sim_data1[,2] <= LSC & sim_data1[,2] >= -LSC))
      z1=z1+controle_prev1
    }
    
    controle_prevv1 <-(sim_data1[,3] <= LSC & sim_data1[,3] >= -LSC)
    
    if (controle_prevv1==0){
      D=c(D,1)
      s1=s1+1
      GG[i]=sum(D)
      next
    }
    
    D=c(D,1)
    
    z2=0
    while (z2==0){
      sim_data2 <- rmvnorm(1, mean = mu3b, sigma = Sigma3)
      controle_prev2 <- ((sim_data2[,1] <= LSC & sim_data2[,1] >= -LSC) &
                           (sim_data2[,2] <= LSC & sim_data2[,2] >= -LSC))
      z2=z2+controle_prev2
    }
    
    controle_prevv2 <-(sim_data2[,3] <= LSC & sim_data2[,3] >= -LSC)
    
    if (controle_prevv2==0){
      D=c(D,1)
      s1=s1+1
      GG[i]=sum(D)
      next
    }
    D=c(D,1)
    
    while (s1<1){
      
      z3=0
      while (z3==0){
        sim_data3 <- rmvnorm(1, mean = mu3c, sigma = Sigma3)
        controle_prev3 <- ((sim_data3[,1] <= LSC & sim_data3[,1] >= -LSC) &
                             (sim_data3[,2] <= LSC & sim_data3[,2] >= -LSC))
        z3=z3+controle_prev3
      }
      controle_prevv3 <-(sim_data3[,3] <= LSC & sim_data3[,3] >= -LSC)
      
      if (controle_prevv3==0){
        D=c(D,1)
        s1=s1+1
      }else{
        D=c(D,1)
        
      }
      
    }
    
    GG[i]=sum(D)
  }
  ARL1=mean(GG)
  SDRL=(var(GG))^0.5
  Result[ii,1]=u1
  Result[ii,2]=ARL1
  Result[ii,3]=SDRL
}
Result=data.frame(Result)
writexl::write_xlsx(Result,"Result_Sim.xlsx")
toc()