library(iterators)
library(foreach)
library(doParallel)

p=10 # number of dimension
n=100 # number of observations 
M=10  # number of machine
c=1 # the catagory of precision matrix
nworkers <- detectCores()-1 # number of cores

X=data_generate(p,n/2,M,c)
Y=data_generate(p,n/2,M,c)
Omega=X$matrix_Omega;
local_locations_aggre1=list();
local_locations_aggre2=list();
cl <- makePSOCKcluster(nworkers)
clusterSetRNGStream(cl, c(1:nworkers))
registerDoParallel(cl)
local_locations_aggre1 = foreach(i = 1:M, .packages = c("Matrix", "glmnet", "MASS")) %dopar%{
  X_temp=matrix(0,n/2,p)
  X_temp=X$Matrix[(1+n*(i-1)/2):(n*i/2), ]
  local_nonzero_locations(X_temp)
}
local_locations_aggre2 = foreach(i = 1:M, .packages = c("Matrix", "glmnet", "MASS")) %dopar%{
  Y_temp=matrix(0,n/2,p)
  Y_temp=Y$Matrix[(1+n*(i-1)/2):(n*i/2), ]
  local_nonzero_locations(Y_temp)
}
stopCluster(cl)
global_locations1=list();
global_locations2=list(); 
global_locations1=center_nonzero_locations(local_locations_aggre1);
global_locations1=global_locations1|t(global_locations1);
Beta1 = matrix(0, nrow = sum(global_locations1), ncol = M )
Rij1 = matrix(0, nrow = sum(global_locations1), ncol = M)
Rii1 = matrix(0, nrow = p, ncol = M)
loc1=as.matrix(which(global_locations1!=0,arr.ind=TRUE))
#loc = as.matrix(summary(global_locations)[1:2])
rbeta1= matrix(0, nrow = sum(global_locations1), ncol = M )
for (i in 1:M){
  temp1=local_est(Y$Matrix[((i-1)*n/2+1):(n*i/2),],global_locations1);
  rbeta1[,i] = summary(t(temp1$local_beta) * temp1$local_cor_diag  +   t(t(temp1$local_beta) * temp1$local_cor_diag))[,3]
  Beta1[,i] =summary(temp1$local_beta)[,3]
  Rii1[,i] = temp1$local_cor_diag
  Rij1[,i] =summary(temp1$local_cor)[,3]
}
Omega_est1=global_est(loc1,Blocal = Beta1, rijlocal = Rij1,riilocal = Rii1,rbeta =  rbeta1, global_locations1)
global_locations2=center_nonzero_locations(local_locations_aggre2);
global_locations2=global_locations2|t(global_locations2);
Beta2 = matrix(0, nrow = sum(global_locations2), ncol = M )
Rij2 = matrix(0, nrow = sum(global_locations2), ncol = M)
Rii2 = matrix(0, nrow = p, ncol = M)
loc2=as.matrix(which(global_locations2!=0,arr.ind=TRUE))
#loc = as.matrix(summary(global_locations)[1:2])
rbeta2= matrix(0, nrow = sum(global_locations2), ncol = M )
for (i in 1:M){
  temp2=local_est(X$Matrix[((i-1)*n/2+1):(n*i/2),],global_locations2);
  rbeta2[,i] = summary(t(temp2$local_beta) * temp2$local_cor_diag  +   t(t(temp2$local_beta) * temp2$local_cor_diag))[,3]
  Beta2[,i] =summary(temp2$local_beta)[,3]
  Rii2[,i] = temp2$local_cor_diag
  Rij2[,i] =summary(temp2$local_cor)[,3]
}
Omega_est2=global_est(loc2,Blocal = Beta2, rijlocal = Rij2,riilocal = Rii2,rbeta =  rbeta2, global_locations2)
Omega_est=(Omega_est1+Omega_est2)/2
print(Omega)
print(Omega_est)