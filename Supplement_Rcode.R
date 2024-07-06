library(Rcpp)
library(matlib)
library(reshape2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
path = getwd()
sourceCpp(paste(path,"/CAR_cpp.cpp",sep=''))
# Note that the sourceCPP file include codes of STRPB, PS, HH randomization procedures

expit = function(x){
  y = 1/(1+exp(-x))
  return(y)
}

multi_level_to_single_level<-function(z){
  n_level<-dim(z)[1]
  n<-dim(z)[2]
  fz<-numeric(n)
  for(i in 1:n_level){
    fz<-fz+z[i,]*i
  }
  return(fz)
}

single_level_to_multi_level<-function(z){
  n_level<-length(unique(z))
  n = length(z)
  fz<-matrix( numeric(n_level*n), nrow = n_level)
  for(i in 1:n_level){
    fz[i,]<-as.numeric(z == i)
  }
  return(fz)
}


########## Complete randomization ##########

CR   <- function( n , rho = 0.5 ){
  Treatment    <- rbinom( n , 1 ,rho)
  return( Treatment )
}


########### Treatment assignment, CR, PS, HH, STRPB ##############

treatment_assignment<-function(n,X,randomization,rho = 1/2){
  wt = rep(1,ncol(X)) # weight of margin in PS
  p=0.85
  wt2 = c(0.1,rep(0.4/ncol(X),ncol(X)),0.5) # weight  in HH
  k = 2
  BS = 4 # Block size
  if(randomization=="CR"){
    I<-CR(n, rho)
  }else if(randomization=="PS"){
    if(rho!=1/2){stop("For now, the proportion of treatment under CABC has to be 1/2.")}
    I<-PS_cpp(X,wt,rho,p)
  }else if(randomization=="STRPB"){
    if(rho!=1/2){stop("For now, the proportion of treatment under permuted_block has to be 1/2.")}
    I<-STRPB_cpp(X,k,BS,rho)
  }else if(randomization=="HH"){
    if(rho!=1/2){stop("For now, the proportion of treatment under CABC has to be 1/2.")}
    I<-HH_cpp(X,wt2,rho,p)
  }
  return(I)
}


########### First simulation data based on four randomization procedure ##############


data_gen<-function(n,theta,randomization = c('CR', 'PS', 'STRPB', 'HH'), rho,case=c("case1")){
  ################ Case 1, no time-variant coefficient, full model ##################
  if(case=="case1"){ # multiple z, function of t is linear
    beta0 = -1
    beta1 = 1
    beta2 = 1.5
    beta3 = 2
    z1<-rmultinom(n,1,c(0.5,0.5))
    z2<-rmultinom(n,1,c(0.4,0.3,0.3))
    Z1 = multi_level_to_single_level(z1)
    Z2 = multi_level_to_single_level(z2)
    X = data.frame(Z1,Z2)
    X = apply(X, 2, as.numeric) 
    n_strata<-dim(z1)[1]*dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(z1)[1]){
      for(j in 1:dim(z2)[1]){
        strata_z[,ind]<-z1[i,]*z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,X,randomization,rho)
    
    eta = beta0 + theta*I + beta1*z1[1,]+ beta2*z2[1,]+ beta3*z2[2,] 
    p = 1/(1+exp(-eta))
    y = rbinom(n,1,p)
    data.simu<-data.frame(y,I, 1-I,z1[1,],z2[1,],z2[2,],z1[1,],z2[1,],z2[2,],strata_z)
    names(data.simu)<- c("y", "I1", "I0","model_Z1","model_Z21","model_Z22","model_R1","model_R21","model_R22",paste0("strata",1:n_strata))
    return(data.simu)
  }
  if(case=="case2"){ # multiple z, function of t is linear
    beta0 = -1
    beta1 = 1
    beta2 = 1.5
    beta3 = 2
    z1 = rnorm(n,0,1)
    z2 = rnorm(n,0,1)
    z1_cat<- 2 - as.numeric(z1< qnorm(0.5))
    z2_cat = numeric(n)
    z2_cat[z2 < qnorm(0.4)] = 1
    z2_cat[qnorm(0.4) <= z2 &  z2< qnorm(0.7)] = 2
    z2_cat[z2 >= qnorm(0.7)] = 3
    
    Z1 = single_level_to_multi_level(z1_cat)
    Z2 = single_level_to_multi_level(z2_cat)
    
    X = data.frame(z1_cat,z2_cat)
    X = apply(X, 2, as.numeric) 
    n_strata<-dim(Z1)[1]*dim(Z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(Z1)[1]){
      for(j in 1:dim(Z2)[1]){
        strata_z[,ind]<-Z1[i,]*Z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,X,randomization,rho)
    
    eta = beta0 + theta*I + beta1*z1+ beta2*z2 
    p = 1/(1+exp(-eta))
    y = rbinom(n,1,p)
    data.simu<-data.frame(y,I, 1-I,z1,z2,z1,z2,strata_z)
    names(data.simu)<- c("y", "I1", "I0","model_Z1","model_Z2","model_R1","model_R2",paste0("strata",1:n_strata))
    return(data.simu)
  }
  if(case=="case3"){ # multiple z, function of t is linear
    z1<-rmultinom(n,1,c(0.5,0.5))
    z2<-rmultinom(n,1,c(0.5,0.5))
    z3<-rmultinom(n,1,c(0.5,0.5))
    z4<-rmultinom(n,1,c(0.5,0.5))
    z5<-rmultinom(n,1,c(0.5,0.5))
    z6<-rmultinom(n,1,c(0.5,0.5))
    z7<-rmultinom(n,1,c(0.5,0.5))
    
    
    Z1 = multi_level_to_single_level(z1)
    Z2 = multi_level_to_single_level(z2)
    Z3 = multi_level_to_single_level(z3)
    Z4 = multi_level_to_single_level(z4)
    Z5 = multi_level_to_single_level(z5)
    Z6 = multi_level_to_single_level(z6)
    Z7 = multi_level_to_single_level(z7)
    
    X = data.frame(Z1,Z2,Z3,Z4,Z5,Z6,Z7)
    X = apply(X, 2, as.numeric) 
    n_strata<-dim(z1)[1]*dim(z2)[1]*dim(z3)[1]*dim(z4)[1]*dim(z5)[1]*dim(z6)[1]*dim(z7)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i1 in 1:dim(z1)[1]){
      for(i2 in 1:dim(z2)[1]){
        for(i3 in 1:dim(z3)[1]){
          for(i4 in 1:dim(z4)[1]){
            for(i5 in 1:dim(z5)[1]){
              for(i6 in 1:dim(z6)[1]){
                for(i7 in 1:dim(z7)[1]){
                  strata_z[,ind]<-z1[i1,]*z2[i2,]*z3[i3,]*z4[i4,]*z5[i5,]*z6[i6,]*z7[i7,]
                  ind<-ind+1
                }
              }
            }
          }
        }
      }
    }
    I<-treatment_assignment(n,X,randomization,rho)
    X = X - 1
    beta0 = -3
    beta = rep(1,7)
    eta = beta0 + theta*I + + as.vector(beta%*%t(X))
    p = 1/(1+exp(-eta))
    y = rbinom(n,1,p)
    ### multiplicative Z, first combination means covariate in working model, next combination means covaraite in randomization
    data.simu<-data.frame(y,I, 1-I,X,X, strata_z)
    names(data.simu)<- c("y", "I1", "I0",paste0("model_Z",1:7),paste0("model_R",1:7),paste0("strata",1:n_strata))
    return(data.simu)
  }
}

########### Second simulation data based on four randomization procedure ##############

data_gen2<-function(n,theta,randomization = c('CR', 'PS', 'STRPB', 'HH'), rho,case=c("case1")){
  ################ Case 1, no time-variant coefficient, full model ##################
  if(case=="case1"){ # multiple z, function of t is linear
    beta0 = -1
    beta1 = 1
    beta2 = 1
    z1<-rmultinom(n,1,c(0.5,0.5))
    z2<-rmultinom(n,1,c(0.4,0.3,0.3))
    Z1 = multi_level_to_single_level(z1)
    Z2 = multi_level_to_single_level(z2)
    X = data.frame(Z1,Z2)
    X = apply(X, 2, as.numeric) 
    n_strata<-dim(z1)[1]*dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(z1)[1]){
      for(j in 1:dim(z2)[1]){
        strata_z[,ind]<-z1[i,]*z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,X,randomization,rho)
    Z1 = scale(Z1, scale = FALSE)
    Z2 = scale(Z2, scale = FALSE)
    
    eta = beta0 + theta*I + beta1*exp(Z1)+ beta2*exp(Z2)
    p = 1/(1+exp(-eta))
    y = rbinom(n,1,p)
    data.simu<-data.frame(y,I, 1-I,z1[1,],z2[1,],z2[2,],strata_z)
    names(data.simu)<- c("y", "I1", "I0","model_Z1","model_Z21","model_Z22",paste0("strata",1:n_strata))
    return(data.simu)
  }
  
  if(case=="case2"){ # multiple z, function of t is linear
    beta0 = -1
    beta1 = 1
    beta2 = 1
    z1 = rnorm(n,0,1)
    z2 = rnorm(n,0,1)
    z1_cat<- 2 - as.numeric(z1< qnorm(0.5))
    z2_cat = numeric(n)
    z2_cat[z2 < qnorm(0.4)] = 1
    z2_cat[qnorm(0.4) <= z2 &  z2< qnorm(0.7)] = 2
    z2_cat[z2 >= qnorm(0.7)] = 3
    
    Z1 = single_level_to_multi_level(z1_cat)
    Z2 = single_level_to_multi_level(z2_cat)
    
    X = data.frame(z1_cat,z2_cat)
    X = apply(X, 2, as.numeric) 
    n_strata<-dim(Z1)[1]*dim(Z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(Z1)[1]){
      for(j in 1:dim(Z2)[1]){
        strata_z[,ind]<-Z1[i,]*Z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,X,randomization,rho)
    
    eta = beta0 + theta*I + beta1*exp(z1)+ beta2*z2^2
    p = 1/(1+exp(-eta))
    y = rbinom(n,1,p)
    data.simu<-data.frame(y,I, 1-I,strata_z)
    names(data.simu)<- c("y", "I1", "I0",paste0("strata",1:n_strata))
    return(data.simu)
  }
  if(case=="case3"){ # multiple z, function of t is linear
    beta0 = -1
    beta1 = 1
    beta2 = 1
    beta3 = 1
    z1 = rnorm(n,0,1)
    z2 = rnorm(n,0,1)
    z1_cat<- 2 - as.numeric(z1< qnorm(0.5))
    z2_cat = numeric(n)
    z2_cat[z2 < qnorm(0.4)] = 1
    z2_cat[qnorm(0.4) <= z2 &  z2< qnorm(0.7)] = 2
    z2_cat[z2 >= qnorm(0.7)] = 3
    
    Z1 = single_level_to_multi_level(z1_cat)
    Z2 = single_level_to_multi_level(z2_cat)
    
    X = data.frame(z1_cat,z2_cat)
    X = apply(X, 2, as.numeric) 
    n_strata<-dim(Z1)[1]*dim(Z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(Z1)[1]){
      for(j in 1:dim(Z2)[1]){
        strata_z[,ind]<-Z1[i,]*Z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,X,randomization,rho)
    
    eta = beta0 + theta*I + beta1*exp(z1+z2)+ beta2*z2^2 + beta3* log(abs(z1))
    p = 1/(1+exp(-eta))
    y = rbinom(n,1,p)
    data.simu<-data.frame(y,I, 1-I,strata_z)
    names(data.simu)<- c("y", "I1", "I0",paste0("strata",1:n_strata))
    return(data.simu)
  }
}





###############. True marginal RR based on parameter setting ###########################
RR_true = function(n,theta,case = 'case1'){
if(case == 'case1'){
  beta0 = -1
  beta1 = 1
  beta2 = 1.5
  beta3 = 2
  p1 = (0.5*0.4 * expit(beta0 + theta + beta2) + 
          0.5*0.3 * expit(beta0 + theta) + 
          0.5*0.3 * expit(beta0 + theta + beta3) + 
          0.5*0.4 * expit(beta0 + theta + beta1 + beta2) + 
          0.5*0.3 * expit(beta0 + theta + beta1) + 
          0.5*0.3 * expit(beta0 + theta + beta1+ beta3)  )
  
  p0 = (0.5*0.4 * expit(beta0 + beta2) + 
          0.5*0.3 * expit(beta0) + 
          0.5*0.3 * expit(beta0 + beta3) + 
          0.5*0.4 * expit(beta0 + beta1 + beta2) + 
          0.5*0.3 * expit(beta0 + beta1) + 
          0.5*0.3 * expit(beta0 + beta1+ beta3)  )
  True_RR = p1/p0
  return(True_RR )
}
}




########## Traditional Wald test T_L ##################

LogRR_regular = function(logRR, data.simu){
  y1 = data.simu$y[data.simu$I1 == 1]
  y0 = data.simu$y[data.simu$I1 == 0]
  n1 = sum(data.simu$I1 == 1)
  n0 = sum(data.simu$I0 == 1)
  t_test_omit =(log(mean(y1)/mean(y0)) - logRR)/sqrt((1-mean(y1))/mean(y1)/n1 + (1-mean(y0))/mean(y0)/n0)
  return(t_test_omit)
}

########## Adjusted Wald test T_{RL} using standardized variance estimator  ##################

LogRR_adj_unpool = function(logRR,data.simu,randomization){
  y = data.simu$y
  I1 = data.simu$I1
  data_model<-data.simu[,grepl("model_Z",names(data.simu))] 
  dat = data.frame(y, I1, data_model)
  model = glm(y ~  . , data = dat, family = binomial ( link = logit ))
  f = summary(model)
  
  theta_est = f$coefficients[2,1]
  beta0 = f$coefficients[1,1]
  beta = f$coefficients[-c(1,2),1]
  
  S = f$cov.unscaled
  if(identical(beta, numeric(0))){
    log_RR_info = RR_no_covariates(beta0,theta_est,S)
  } else{
    log_RR_info = RR_marginal(data_model,beta0,beta,theta_est,S)
  }
  log_RR_est = log_RR_info[1]
  
  n = length(y)
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  n_strata<-dim(strata_z)[2]
  z<-strata_z
  var_p1<-numeric(n_strata)
  var_p0<-numeric(n_strata)
  E_p1<-numeric(n_strata)
  E_p0<-numeric(n_strata)
  prob_z<-numeric(n_strata)
  n_s<-numeric(n_strata)
  
  for(j in 1:n_strata){
    var_p1[j]<-var(y[z[,j]==1 & I1==1])
    var_p0[j]<-var(y[z[,j]==1 & I1==0])
    E_p1[j]<-mean(y[z[,j]==1 & I1==1])
    E_p0[j]<-mean(y[z[,j]==1 & I1==0])
    prob_z[j]<-mean(z[,j])
    n_s[j] = sum(z[,j] ==1)
  }
  var_p1[is.na(var_p1)] = 0
  var_p0[is.na(var_p0)] = 0
  ##########
  var_H1 = (var_p1)%*%prob_z/(mean(y[I1==1])^2) + (var_p0)%*%prob_z/(mean(y[I1==0])^2)
  ########
  var_H2 = ((E_p1 - mean(y[I1==1]))/mean(y[I1==1]) - (E_p0 - mean(y[I1==0]))/mean(y[I1==0]) )^2%*%prob_z
  var_H3 = ((E_p1 - mean(y[I1==1]))/mean(y[I1==1]) + (E_p0 - mean(y[I1==0]))/mean(y[I1==0]) )^2%*%prob_z
  
  if(randomization == 'CR'){
    var_adj = 2*var_H1+ var_H2 + var_H3
  }
  else{
    var_adj = 2*var_H1+ var_H2
  }
  
  
  t_test_adj =   (log_RR_est - logRR)/(sqrt(var_adj)/sqrt(n))
  return(t_test_adj)
}


############# Estimated marginal log(RR) and its variance adjusted on logistic regression model##############
##### x -- covariate vector, beta_0, beta, theta -- estimated parameters based on working logistic model,
####### S -- estimated variance based on working logistic model #######################


RR_marginal = function(x,beta0,beta,theta,S){
  x = as.matrix(x)
  k = dim(x)[2]
  p1 = 1/(1 + exp( -(theta + beta%*%t(x) + beta0) ) ) 
  p0 = 1/(1 + exp( -( beta%*%t(x) + beta0) ) )
  p1.x = matrix(p1*(1-p1)*x[,1], nrow = nrow(x))
  p0.x = matrix(p0*(1-p0)*x[,1], nrow = nrow(x))
  if(k>=2){
    for(i in 2:k){
      p1.x = cbind(p1.x,t(p1*(1-p1)*x[,i]))
      p0.x = cbind(p0.x, t(p0*(1-p0)*x[,i]))
    }
  }
  
  log_RR_margin = log(mean(p1)/mean(p0))
  
  gradient_f = 1/mean(p1) * c( mean(p1*(1-p1)), mean(p1*(1-p1)), colMeans(p1.x) ) - 
    1/mean(p0) * c( mean(p0*(1-p0)),0, colMeans(p0.x) )
  Var_log_RR_margin =  gradient_f %*% S %*% gradient_f
  
  return(c(log_RR_margin,Var_log_RR_margin) )
}

########## Estimated marginal log(RR) when no covariates are included in the working model  ############
##### beta_0,theta -- estimated parameters based on working logistic model,
####### S -- estimated variance based on working logistic model #######################

RR_no_covariates = function(beta0,theta,S){
  p1 = 1/(1 + exp( -(theta + beta0) ) ) 
  p0 = 1/(1 + exp( -(  beta0) ) )
  log_RR_margin = log(mean(p1)/mean(p0))
  
  gradient_f = 1/mean(p1) * c( mean(p1*(1-p1)), mean(p1*(1-p1)) ) - 
    1/mean(p0) * c(  mean(p0*(1-p0)),0)
  Var_log_RR_margin =  gradient_f %*% S %*% gradient_f
  
  return(c(log_RR_margin,Var_log_RR_margin) )
}



########## Traditional Wald test $T_S$ based on logistic regression ##################
######## logRR_0: null hypothesis value log(RR_0), set as 0 in Type-I error test example #########


logistic_test = function(logRR_0, data.simu){
  y = data.simu$y
  I1 = data.simu$I1
  data_model<-data.simu[,grepl("model_Z",names(data.simu))] 
  dat = data.frame(y, I1, data_model)
  model = glm(y ~  . , data = dat, family = binomial ( link = logit ))
  f = summary(model)
  
  theta_est = f$coefficients[2,1]
  beta0 = f$coefficients[1,1]
  beta = f$coefficients[-c(1,2),1]
  
  S = f$cov.unscaled
  if(identical(beta, numeric(0))){
    log_RR_info = RR_no_covariates(beta0,theta_est,S)
  } else{
    log_RR_info = RR_marginal(data_model,beta0,beta,theta_est,S)
  }
  log_RR_est = log_RR_info[1]
  log_RR_sd = sqrt(log_RR_info[2])
  
  t_test = (log_RR_est - logRR_0)/log_RR_sd
  
  return(t_test)
}



################### adjusted test test $T_{RS}$ using standardization  variance estimator ###############

LogRR_adj = function(logRR,data.simu, randomization){
  y = data.simu$y
  I1 = data.simu$I1
  data_model<-data.simu[,grepl("model_Z",names(data.simu))] 
  dat = data.frame(y, I1, data_model)
  model = glm(y ~  . , data = dat, family = binomial ( link = logit ))
  f = summary(model)
  
  theta_est = f$coefficients[2,1]
  beta0 = f$coefficients[1,1]
  beta = f$coefficients[-c(1,2),1]
  S= f$cov.unscaled
  
  Theta = c(beta0, theta_est, beta)
  
  log_RR_info = RR_marginal(data_model,beta0,beta,theta_est,S)
  log_RR_est = log_RR_info[1]
  
  
  X = dat
  X[,1] =  rep(1,n)
  
  h_x = expit(as.matrix(X)%*%as.matrix(Theta))
  G = (y - h_x)*X
  
  n_col = length(Theta)
  
  
  p1 = h_x *(1-h_x)*dat[,3:n_col]
  Theta_0 = c(beta0, 0, beta)
  h0_x = expit(as.matrix(X)%*%as.matrix(Theta_0))
  p0 =  h0_x *(1-h0_x)*X[,3:n_col]
  gradient_f = 1/mean(h_x) * c(mean(h_x*(1-h_x)), mean(h_x*(1-h_x)), colMeans(p1) ) - 
    1/mean(h0_x) * c(  mean(h0_x*(1-h0_x)),0, colMeans(p0) )
  O = as.vector(gradient_f %*%   S %*% t(G))*sqrt(n)
  O1 = O*I1
  O0 = O*(1-I1)

  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  n_strata<-dim(strata_z)[2]
  z<-strata_z
  var_p<-numeric(n_strata)
  var_p1<-numeric(n_strata)
  var_p0<-numeric(n_strata)
  
  E_p1 = numeric(n_strata)
  E_p0 = numeric(n_strata)
  
  prob_z<-numeric(n_strata)
  n_s<-numeric(n_strata)
  for(j in 1:n_strata){
    var_p[j]<-var(O[z[,j]==1])
    var_p1[j]<-var(O1[z[,j]==1 & I1==1])
    var_p0[j]<-var(O[z[,j]==1 & I1==0])
    E_p1[j]<-mean(O[z[,j]==1 & I1==1])
    E_p0[j]<-mean(O[z[,j]==1 & I1==0])
    prob_z[j]<-mean(z[,j])
    n_s[j] = sum(z[,j] ==1)
  }
  var_p1[is.na(var_p1)] = 0
  var_p0[is.na(var_p0)] = 0
  
  Sigma1 = 1/2*(var_p1 +var_p0 )%*%prob_z
  Sigma2 =  1/4*((E_p1 - mean(O[I1==1])) + (E_p0 - mean(O[I1==0]) ))^2%*%prob_z 
  Sigma3 =  1/4*((E_p1 - mean(O[I1==1])) - (E_p0 - mean(O[I1==0]) ))^2%*%prob_z 
  
  if(randomization == 'CR'){
    var_adj =   Sigma1 +  Sigma2 + Sigma3 
  }
  else{
    var_adj =  Sigma1 +  Sigma2
  }
  se_adj = sqrt(var_adj)
  t_test_adj =   (log_RR_est - logRR)/ se_adj
  return(t_test_adj)
}

################### adjusted test test T^*_{RS} and T^*_{RL} using  variance estimator from full logistic model  ###############

LogRR_adj_parametric = function(logRR,data.simu,randomization){
  y = data.simu$y
  I1 = data.simu$I1
  data_model<-data.simu[,grepl("model_Z",names(data.simu))] 
  dat = data.frame(y, I1, data_model)
  model = glm(y ~  . , data = dat, family = binomial ( link = logit ))
  f = summary(model)
  
  theta_est = f$coefficients[2,1]
  beta0 = f$coefficients[1,1]
  beta = f$coefficients[-c(1,2),1]
  
  S = f$cov.unscaled
  if(identical(beta, numeric(0))){
    log_RR_info = RR_no_covariates(beta0,theta_est,S)
  } else{
    log_RR_info = RR_marginal(data_model,beta0,beta,theta_est,S)
  }
  log_RR_est = log_RR_info[1]
  if(randomization == 'CR'){
    log_RR_sd = sqrt(log_RR_info[2])
    t_test_adj =   (log_RR_est - logRR)/log_RR_sd
  }
  else{
  data_model_F<-data.simu[,grepl("model_R",names(data.simu))] 
  dat_F = data.frame(y, I1, data_model_F)
  model_F = glm(y ~  . , data = dat_F, family = binomial ( link = logit ))
  f_F = summary(model_F)
  
  theta_est = f_F$coefficients[2,1]
  beta0 = f_F$coefficients[1,1]
  beta = f_F$coefficients[-c(1,2),1]
  
  S = f_F$cov.unscaled
  if(identical(beta, numeric(0))){
    log_RR_info = RR_no_covariates(beta0,theta_est,S)
  } else{
    log_RR_info = RR_marginal(data_model_F,beta0,beta,theta_est,S)
  }
  log_RR_sd = sqrt(log_RR_info[2])
  
  t_test_adj =   (log_RR_est - logRR)/log_RR_sd
  }
  return(t_test_adj)
}



################################### exam test on type-I error and power ########
################# return value: T_F, T_P, T_0 means wald test on full model, partial model and zero-covariate model; 
############### T_P_adj, T_0_adj means adjusted wald test $T_{RS}, T_{RL} in partial and zero-covariate model ##################
############### T_P_adj_parm, T_0_adj_parm means adjusted wald test $T^*_{RS}, T^*_{RL} in partial and zero-covariate model ##################



exam_power<-function(n, theta, logRR = 0, randomization, rho,case = 'case1',alpha = 0.05,n_rep = 10){
  rej_count<-numeric(8)
  for(k in 1:n_rep){
    print(k)
    data.simu<-data_gen(n, theta, randomization, rho,case )
    if(case == 'case1'){
      data.simu_p = data.simu[,c(-5,-6)]
      data.simu_0 = data.simu[,c(-4,-5,-6)]
    }
    else if(case == 'case2'){
      data.simu_p = data.simu[,c(-5)]
      data.simu_0 = data.simu[,c(-4,-5)]
    }
    else if(case == 'case3'){
      data.simu_p = data.simu[,c(-7,-8,-9,-10)]
      data.simu_0 = data.simu[,c(-4,-5,-6,-7,-8,-9,-10)]
    }
    ########## Full model #############
    t_f =  logistic_test(logRR, data.simu)
    t_f_adj =  LogRR_adj(logRR, data.simu,randomization)
    if(abs(t_f)>abs(qnorm(alpha/2))){
      rej_count[1]<-rej_count[1]+1
    }
    if(abs(t_f_adj)>abs(qnorm(alpha/2))){
      rej_count[2]<-rej_count[2]+1
    }
    ################ Partial Model #############    
    t_p<-logistic_test(logRR, data.simu_p)
    t_p_adj =  LogRR_adj(logRR, data.simu_p, randomization)
    t_p_adj_parm = LogRR_adj_parametric(logRR, data.simu_p,randomization)
    
    if(abs(t_p)>abs(qnorm(alpha/2))){
      rej_count[3]<-rej_count[3]+1
    }
    if(abs(t_p_adj)>abs(qnorm(alpha/2))){
      rej_count[4]<-rej_count[4]+1
    }
    if(abs(t_p_adj_parm)>abs(qnorm(alpha/2))){
      rej_count[5]<-rej_count[5]+1
    }
    
    ################ Zero-covariate Model ############# 
    
    ####### No adjustment #######
    t_0<-LogRR_regular(logRR,data.simu_0)
    
    t_0_adj =  LogRR_adj_unpool(logRR, data.simu_0,randomization)
    t_0_adj_parm = LogRR_adj_parametric(logRR, data.simu_0,randomization)
    
    if(abs(t_0)>abs(qnorm(alpha/2))){
      rej_count[6]<-rej_count[6]+1
    }
    if(abs(t_0_adj)>abs(qnorm(alpha/2))){
      rej_count[7]<-rej_count[7]+1
    }
    if(abs(t_0_adj_parm)>abs(qnorm(alpha/2))){
      rej_count[8]<-rej_count[8]+1
    }
    # 
    print(rej_count)
  }
  power_res<-rej_count/n_rep
  power_res_mat<-matrix(power_res,nrow = 1,byrow = TRUE)
  colnames(power_res_mat)<-c("T_F","T_F_adj","T_P","T_P_adj","T_P_adj_parm", "T_0","T_0_adj","T_0_adj_parm")
  rownames(power_res_mat)<-"empirical_power"
  return(power_res_mat)
}

########## Type I error test ##########

typeI_error = matrix(numeric(4*8), ncol = 8)

n = 500
case = 'case1'
rho = 0.5 
iter = 2000
theta = 0
logRR = 0

set.seed(1234)

#########
case = 'case1'
randomization = 'CR'
typeI_error[1,] = exam_power(n, theta, logRR, randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
randomization = 'STRPB'
typeI_error[2,] = exam_power(n, theta, logRR, randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
randomization = 'PS'
typeI_error[3,] = exam_power(n, theta, logRR, randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
randomization = 'HH'
typeI_error[4,] = exam_power(n, theta, logRR, randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)


################################################
############## Power test ############

set.seed(1234)
theta = seq(0,0.9,0.1)
iter = 2000
test_result_CR = matrix(numeric(10*8), nrow = 10)
test_result_STRPB = matrix(numeric(10*8), nrow = 10)
test_result_PS = matrix(numeric(10*8), nrow = 10)
test_result_HH = matrix(numeric(10*8), nrow = 10)

randomization = 'STRPB'
for(i in 1:length(theta)){
  test_result_STRPB[i,] = exam_power(n, theta[i],logRR, randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
}
randomization = 'PS'
for(i in 1:length(theta)){
  test_result_PS[i,] = exam_power(n, theta[i], logRR,randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
}
randomization = 'HH'
for(i in 1:length(theta)){
  test_result_HH[i,] = exam_power(n, theta[i], logRR,randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
}
randomization = 'CR'
for(i in 1:length(theta)){
  test_result_CR[i,] = exam_power(n, theta[i], logRR,randomization, rho = 1/2,case,alpha = 0.05,n_rep = iter)
}

