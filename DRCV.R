library(MASS)
library(glmnet)
library(Matrix)
library(HDeconometrics)
library(psych)
data_generate <- function(p_value,n_value,M_value,c_value=1){
  #p means number of variables
  #n means number of samples
  #M means number of machine
  #c means category of precision matrix,the value only equals 1, 2, 3
  p<-p_value
  n<-n_value
  M<-M_value
  c<-c_value
  #according the value of c, generate the specific adjacency matrix A
  #when c=1, the graph is banded graph;
  if(c==1){
    matrix_Adja=matrix(0,p,p);
    for(i in 1:p-1){
      for(j in (i+1):p){
        if (abs(i-j)<=3)
        {matrix_Adja[i,j]=1;
        matrix_Adja[j,i]=1;}
        else
        {matrix_Adja[i,j]=0;}
      }
    }
  }
  #when c=2. the graph is Erdos-Renyi random graph;
  if(c==2){
      matrix_Adja=matrix(0,p,p);
      for (i in 1:p-1){
        for (j in (i+1):p){
          matrix_Adja[i,j]=rbinom(1,1,5/p)
          matrix_Adja[j,i]=matrix_Adja[i,j]
        }
      }
      #matrix_Adja=matrix_Adja_case4
  }  
  #when c=3, the graph is cluster graph;
  if(c==3){
      matrix_Adja=matrix(0,p,p);
      count=ceiling(p/10);
      for (d in 1:count){
        if (d!=count){
          for(i in ((d-1)*10+1):(d*10-1)){
            for (j in (i+1):(d*10)){
             matrix_Adja[i,j]=rbinom(1,1,0.6);
              matrix_Adja[j,i]=matrix_Adja[i,j];
            }
          }
        }
        if (d==count){
          for(i in ((d-1)*10+1):(p-1)){
            for (j in (i+1):p){
              matrix_Adja[i,j]=rbinom(1,1,0.6);
              matrix_Adja[j,i]=matrix_Adja[i,j];
            }
          }
        }
    }
    #matrix_Adja=matrix_Adja_case7
  }
  #we reset the non-zero off diagonal elements value to be 0.3;
  matrix_Adja=0.3*matrix_Adja;
  #get the smallest eigenvalue of the adjacency matrix A
  lamada_min=min(eigen(matrix_Adja)$values);
  #define a diagonal matrix D
  matrix_D=matrix(0,p,p);
  for (i in 1:p){
    if (i<p/2){
      matrix_D[i,i]=1;
    }
    else{
      matrix_D[i,i]=1.5;
    }
  }
  #get the precision matrix Omega;
  matrix_Omega=matrix_D%*%(matrix_Adja+(abs(lamada_min)+0.2)*diag(1,p,p))%*%matrix_D;
  #get the covariance matrix Sigma;
  matrix_Sigma=matrix(0,p,p);
  matrix_Sigma=solve(matrix_Omega)
  #get the n i.i.d observations from the multivariate Gaussian distribution
  X=mvrnorm(n*M,rep(0,p),matrix_Sigma);
  return (list(Matrix=X,matrix_Omega = matrix_Omega))
}
local_nonzero_locations <- function(data,smin = 5){
  #data means observation sample
  #p means number of variables
  #n means number of samples
  X=data;
  p=ncol(X);
  n=nrow(X);
  #define a null set beta;
  local_nonzero=sparseMatrix(i=1,j=1,x=0,dims = c(p,p));
  for (k in 1:p){
    #extract the independent variable and dependent variable 
    X_independent_variable=X[,-k];
    Y_dependent_variable=X[,k];
    #use lasso method to estimate parameter_beta.
    #plot(fun_lasso1,xvar="lambda");
    #print(fun_lasso1);
    fun_lasso_fit=cv.glmnet(x=X_independent_variable, y=Y_dependent_variable, family="gaussian", nlambda=100, alpha=1, type.measure = "mse", nfolds = 10)
    #plot(fun_lasso_fit);
    #print(fun_lasso_fit);
    #fun_lasso_fit2 = ic.glmnet(x=X_independent_variable,y=Y_dependent_variable,family="gaussian",nlambda=100,alpha=1,crit = "aic")
    fun_lasso_best=glmnet(x=X_independent_variable, y=Y_dependent_variable, family="gaussian", alpha=1, lambda = fun_lasso_fit$lambda.1se)
    #fun_lasso_best=glmnet(x=X_independent_variable,y=Y_dependent_variable,family="gaussian",alpha=1,lambda = fun_lasso_fit2$lambda)
    #plot(fun_lasso_best);
    #print(fun_lasso_best);
    coef(fun_lasso_best)
    #extract the local machine's nonzero locations and save it in beta
    coe<-as.numeric(coef(fun_lasso_best))
    Active_Index<-which(coe!=0)
    Active_Index=Active_Index[Active_Index>1]-1;
    len=length(Active_Index)
    if(len < smin){
      fun_lasso_more=glmnet(x=X_independent_variable,y=Y_dependent_variable,family="gaussian",alpha=1)
      coe_new = as.numeric(coef(fun_lasso_more)[,which.max(fun_lasso_more$df >=smin)])
      Active_Index<-which(coe_new!=0)
      Active_Index=Active_Index[Active_Index>1]-1;
      local_nonzero[k,c(1:p)[-k][Active_Index]]=1
    }else{
      local_nonzero[k,c(1:p)[-k][Active_Index]]=1
    }
  }
  return(local_nonzero)
  #?̦Z?????????o????x?󾹪??D?s?????
}
center_nonzero_locations <- function (total_nonzero_locations,threshold=0.5){
  #local_locations_aggre means all local machine nonzero locations position aggregation, usually is a list with M items;
  #define total nonzero locations
  M <- length(total_nonzero_locations) # Number of machines
  p <- dim(total_nonzero_locations[[1]])[1] # Number of variables
  #get the center nonzero locations;
  center_nonzero=sparseMatrix(i=1,j=1,x=0,dims = c(p,p));
  for (i in 1:M){
    center_nonzero=center_nonzero+total_nonzero_locations[[i]];
  }
  M_thre=M*threshold;
  center_nonzero=center_nonzero>=M_thre
  return(center_nonzero)
}
local_est <- function(X, global_locations){
  #global_locations means global nonzero location
  #data means observation samples
  #define global_nonzero_location
  #define observations
  #p means dimension number and  n means samples number
  p=ncol(X);
  n=nrow(X);
  #define the coefficient beta
  beta=sparseMatrix(i=1,j=1,x=0,dims = c(p,p));
  #define the epsilon
  epsilon=matrix(0,n,p);
  #use the OLS method to estimate beta and epsilon
  for (k in 1:p){
    y_dependent_variables=X[,k];
    #extract the specific columns
    if(sum(global_locations[k,])>0){
      x_independet_variables=matrix(X[,global_locations[k,]], ncol = sum(global_locations[k,]));
      lm_fit=lm.fit(x=cbind(1,x_independet_variables),y = y_dependent_variables); #as.numeric()
      beta[k,global_locations[k,]]=lm_fit$coefficients[-1]
      epsilon[,k]=lm_fit$residuals;
      # lm_fit=lm.fit(x=x_independet_variables,y = y_dependent_variables); #as.numeric()
      # beta[k,global_locations[k,]]=lm_fit$coefficients
      next;
    }
    else {
      epsilon[,k]=y_dependent_variables-colMeans(X)[k]
    }
  }
  #calculate r_{ij}
  r = sparseMatrix(i=1, j=1, x=0, dims = c(p,p));
  rii = rep(0,p);
  for(k in 1:p){
    r[k,global_locations[k,]]=colSums(epsilon[,k]*matrix(epsilon[,global_locations[k,]], ncol = sum(global_locations[k,]) ) ) /n;
    rii[k]=sum(epsilon[,k]*epsilon[,k])/n;
  }  
  result=list(local_beta=beta,local_cor=r,local_cor_diag = rii);
  return(result)
}
global_est <- function(loc, Blocal, rijlocal, riilocal, rbeta, global_locations){
  #local_beta_aggre means all local machine beta aggregation, usually is a 3-dimensional array;  
  #local_r_aggre means all local machine r aggregation, usually is a 3-dimensional array; 
  p <- dim(riilocal)[1];
  M <- dim(riilocal)[2];
  total_beta_aggre = rowMeans(Blocal)
  hat_Tii = rowMeans(riilocal)
  rbetabar = total_beta_aggre * hat_Tii[loc[,2]]
  temp_M = sparseMatrix(i = loc[,1],j = loc[,2], x = rbetabar,dims = c(p,p))
  
  #rbetabar2 = total_beta_aggre * hat_Tii[loc[,2]]
  #temp_M2 = sparseMatrix(i = loc[,1],j = loc[,2], x = rbetabar,dims = c(p,p))
  temp_M = temp_M + t(temp_M)
  
  hat_Tij= rowMeans(rijlocal) - rowMeans(rbeta) +  summary(temp_M)[,3]
  Tpre = sparseMatrix(i = loc[,1],j = loc[,2], x = hat_Tij,dims = c(p,p))
  diag(Tpre) = hat_Tii 
  Omega_pre = Tpre
  Omega_pre = t(t(Omega_pre / hat_Tii)  / hat_Tii)
  return(Omega_pre)
  #define hat_Omega
  #  hat_Omega=sparseMatrix(i=1,j=1,x=0,dims = c(p,p));
  #  for (k in 1:p){
  #    hat_Omega[k,global_locations[k,]]=(hat_Tij[k,global_locations[k,]]/(hat_Tij[k,k])*(1/diag(hat_Tij[global_locations[k,],global_locations[k,]])))
  #    hat_Omega[k,k]=1/hat_Tii[k,k]
  #  }
  #  return(hat_Omega)
}