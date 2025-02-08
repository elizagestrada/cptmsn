
install.packages(c("sn","nortest", "expm"))
library(sn)
library(expm)
library(nortest)

#Estimated canonical form
canonical_form<-function(y){
  n<-NA
  p<-NA
  n<-nrow(y)
  p<-ncol(y)
  
  #Parameter estimation
  param<-msn.mle(y=y)
  xi<-as.vector(param$dp$beta)
  Omega<-param$dp$Omega
  alpha<-param$dp$alpha
  
  omega_inv<-diag(1/sqrt(diag(Omega)))
  Omega_bar<-omega_inv%*%Omega%*%omega_inv
  omega<-diag(sqrt(diag(Omega)))
  delta<-Omega_bar %*% (alpha) / as.numeric(sqrt(1 +  alpha  %*% Omega_bar %*% (alpha)))
  mu_z<-sqrt(2/pi)*delta
  Sigma_z<-Omega_bar-(mu_z)%*%t(mu_z)
  Sigma_y<-omega%*%Sigma_z%*%omega
  
  #Square root matrix of "Inverse Omega"
  sqrt.inv.Omega<-sqrtm(solve(Omega))
  
  #M matrix
  M<-sqrt.inv.Omega%*%Sigma_y%*%sqrt.inv.Omega
  
  #Spectral decomposition of M
  eigen.M<-eigen(M,symmetric = TRUE)
  
  #Q matrix
  Q<-as.matrix(eigen.M$vectors)
  
  #H matrix
  H<-sqrt.inv.Omega%*%Q
  
  
  z <- matrix(rep(NA, p*n), nrow = n, ncol = p)
  for(i in 1:p){
    z[ , i] <- y[ , i] - xi[i]  
  }
  
  #Canonical form
  Z.star <-t( t(H) %*%  t(z))
  return(Z.star)  
  
}

# Test for the univariate skew-normal distribution
Wtest<- function(x){
  n <-length(x)
  cp_w <-sn.mple(y = x, opt.method = "nlminb")$cp
  w <-cp2dp(cp_w, family = "SN")
  vi <-psn(x, w[1], w[2], w[3])   # transformation to uniformity
  vi <-vi[vi != 0]
  vi <-vi[vi != 1]
  yi <-qnorm(vi) # transformation to normality
  sw <-shapiro.test(yi)$statistic # Shapiro-Wilk statistic
  z <-log(n)
  mu_n <- -1.7785544-0.3590050*z -0.072070*z^2 +0.0032656*z^3
  s_n <- 0.1971786 + 0.1653483*z -0.0330235*z^2 + 0.0018982*z^3
  w1 <- log(1-sw)
  pvalue <- pnorm(w1,mu_n, s_n, lower.tail = FALSE)
  return(pvalue)
}


#C test
C_test<-function(y){
  dim<-ncol(y)                #dimension
  can_form<-canonical_form(y) #estimated canonical form
  
  Wpval<-c()                  #p-values of the univariate skewed-normality test
  for (i in 1:dim){
    Wpval[i]<-Wtest(can_form[,i])
  }
  
  #Calculating the Cc-test statistic
  j<-1
  Cc_statistic<-0
  while (j<=dim) {
    Cc_statistic<-Cc_statistic+((1/dim)*tan((0.5-Wpval[j])*pi))  
    j<-j+1
  }
  pval<-pcauchy(Cc_statistic,lower.tail = FALSE)#p-value
  return(pval)  
}
 
#F test
F_test<-function(y){
  dim<-ncol(y)                #dimension
  
  can_form<-canonical_form(y) #estimated canonical form
  
  Wpval<-c()                  #p-values of the univariate skewed-normality test
  for (i in 1:dim){
    Wpval[i]<-Wtest(can_form[,i])
  }
  
  #Calculating the Fc-test statistic
  j<-1
  Fc_statistic<-0
  while (j<=dim) {
    Fc_statistic<-Fc_statistic-2*log(Wpval[j]) 
    j<-j+1
  }
  pval<-pchisq(Fc_statistic,df=2*dim,ncp=0,lower.tail = FALSE)#p-value
  return(pval)
}


#C_d test
Cd_test<-function(y){
  dim<-ncol(y)   #dimension
  
  Wpval<-c()                
  for (i in 1:dim){
    Wpval[i]<-Wtest(y[,i]) #p-values of the univariate skewed-normality test
  }
  
  #Calculating the Cc-test statistic
  j<-1
  Cc_statistic<-0
  while (j<=dim) {
    Cc_statistic<-Cc_statistic+((1/dim)*tan((0.5-Wpval[j])*pi))  
    j<-j+1
  }
  pval<-pcauchy(Cc_statistic,lower.tail = FALSE)#p-value
  return(round(pval,3))
  
}

