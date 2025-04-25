rm(list = ls())
#------------------------------------------------------------------------------------------------------------------
# Simulation specifications (dimensionality):
T = N=n=49 # Number of time points
m = 100 # Number of observation points
m_fore = 110 # Forecast/predict the last m_fore points at state s
N_sim = 100 # Number of simulations

p1=3 ## number of nonzero predictors
p0=1  ## zero predictors
p=p0+p1+1  ##number of covariates,including intercept
corr=0

states = 1:49

# Simulation specifications (parameters):
r_true = 10 # Size parameter for Negative-Binomial; large approximates Poisson
K_true = 4 # True number of factors; 3 or 4 work best
perc_missing = 0.10 # Proportion missing in count observations

W0 = read.csv("D:/49states_adjmatrix.csv")
W = as.matrix(W0[,-1])

# Population:
pop = read.csv("D:/states_population.csv")

# Model specifications:
K = 10 # Number of (unknown) basis functions (K=3 for simulations)



#------------------------------------------------------------------------------------------------------------------
# Load the libraries:
library(dfosr); # devtools::install_github("drkowal/dfosr")
library(coda) # For MCMC output
library(forecast); library(tscount) # For forecasting
library(BayesLogit)
library(MASS)
library(CARBayes)

# Set WD:
setwd("D:/fosr.count_spatial")
source("helper_functions_spatial.R")
source("mcmc_sampler_spatial.R")



# Forecasting study:
#------------------------------------------------------------------------------------------------------------------
# Forecast methods:
names_fore = c('Nb-LCAR', 
               'Pois-LCAR',
               'Guass-LCAR',
               'Nb-IND',
               'Nb-ICAR',
               'Nb-CCAR',
               'Fpc-Basis',
               'Spline-Basis',
               'Spline-NIG',
               'Nb-NIG')  



# Store the errors (compute RMSEs later):
errors_baye =errors_Y =  array(NA, c(length(names_fore), # Number of methods
                                    N_sim))#


# Coverage probability and coverage width:
cover_prob_baye = cover_width_baye = array(NA, c(length(names_fore), N_sim), dimnames = list(names_fore, 1:N_sim))

#############
rmse = matrix(0, length(names_fore), N_sim ) 
alpha_tilde_width =  matrix(0, length(names_fore), N_sim ) 
alpha_tilde_prob =  matrix(0, length(names_fore), N_sim ) 


###############

B_1 = diag(rowSums(W)) - W

V=mvrnorm(1,rep(0,n),0.5^2*ginv(B_1))   ##spatial correlated random component: v_i

L=40   ### L dates for the states
p_star=matrix(0,N,L)
for (i in 1:L) {
  p_star[,i] = exp(mvrnorm(1,V, 0.5^2*diag(N)))
}

p_star= cbind(p_star[,(1:((L/2)-1))] ,0.5*p_star[,((L/2):L)])  ##
prob=p_star/rowSums(p_star)  ## probability


date=NULL   ### date of virus for different states
for(i in 1:N){
  date[i]=which(rmultinom(1, 1,prob[i,])==1)  ### multinomial distribution(L,p)
}


date=readxl::read_excel("D:\\initial date_sim.xlsx")
date = as.matrix(date[,2])
######

for(i in 1:N_sim){
  
  #------------------------------------------------------------------------------------------------------------------
  # Simulate the data:
  #------------------------------------------------------------------------------------------------------------------
  # set.seed(i) # for reproducibility
  sim_data = simulate_fosr.count(m = m, 
                                 W = W, RSNR = 5,
                                 r_true = r_true, 
                                 K_true = K_true, 
                                 p0 = p0, p1 = p1, corr=0,
                                 sparse_factors = FALSE,
                                 perc_missing = perc_missing,
                                 rho=0.7)

  # Store the output:
  Z = sim_data$Z; tau = sim_data$tau ; X=sim_data$X
  
  
  ####  Now delete the observations:
  id_pred=NULL  ### 
  for(ii in 1:n){
    if(date[ii]>(m_fore-m+1))
    { id_pred=rbind(id_pred,ii)
    Z[ii, (m_fore-date[ii]+1):m ]=NA
    }
  }
  
  #  ##### plot
  #   plot(Z[1,],type="l", lwd=3,xlab='Days since the first case', cex.lab=1.3, cex.main=1.5,ylab='Simultaed case counts')
  #   lines(Z[3,],lwd=3,col=3)
  #   lines(Z[4,],lwd=3,col=4)
  #   lines(Z[11,],lwd=3,col=7)
  #lines(Z[12,],lwd=2,col=8)
  #lines(Z[14,],lwd=2,col='chocolate1')
  #   lines(Z[48,],lwd=3,col='chocolate1')
  #lines(Z[49,],lwd=3,col=6)
  #lines(Z[15,],col='chocolate1')
  #lines(Z[43,],lwd=2,col='chocolate1')
  
  #   legend('topleft',ncol=2,legend = c("Alabama",'Arkansas','California','Idaho','Wisconsin'), col=c(1,3,4,7,'chocolate1'),
  #          lwd=3, x.intersp=0.5, text.width=10,seg.len=1)
  
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  # MCMC sampler:  Nb-LCAR
  out_fore_leroux_rho1 = fosr.count_leroux_rho1(Z, tau, X =X, W=W,
                                                K = K, use_approx_poisson = FALSE,
                                                Offset = NULL,
                                                rho.possible = seq(0, 0.95, 0.05),
                                                nsave = 5000, nburn = 5000, nskip = 4,
                                                mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "Pi", "sigma_e","rho")
  )
  
  
  # True mean and true counts (unobserved)
  Z_true=array(NA,c(length(id_pred),m))
  Z_obs =NULL
  
  # Posterior mean and distribution of the forecast:
  Zhat_leroux_rho1 = array(NA,c(length(id_pred),m))
  post_Zpred_leroux_rho1 = NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
     Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
     Z_obs =c(Z_obs,sim_data$Z_full[s,ind_fore])
    
    if(length(ind_fore)>1) 
    {Zhat_leroux_rho1[d,ind_fore] = colMeans(out_fore_leroux_rho1$Zhat[,s,ind_fore])}
    else {Zhat_leroux_rho1[d,ind_fore] = mean(out_fore_leroux_rho1$Zhat[,s,ind_fore]) }
    
    post_Zpred_leroux_rho1=cbind(post_Zpred_leroux_rho1, out_fore_leroux_rho1$Zpred[,s, ind_fore])
  }  
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_leroux_rho1$fk, 
                                              post_alpha_j = out_fore_leroux_rho1$alpha[,j,])
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[1,i] =alpha_tilde_width[1,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[1,i] =alpha_tilde_prob[1,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[1,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[1,i] = alpha_tilde_width[1,i]/((p-1)*m)
  alpha_tilde_prob[1,i] =   alpha_tilde_prob[1,i]/((p-1)*m)
  
  
  ####### mcmc imputation
  #  na.ind = which(is.na(Z), arr.ind = TRUE)
  #  Z_imputed = Z
  #  Z_imputed[na.ind] = rnbinom(n = nrow(na.ind),
  #                              size = out_fore_leroux$r,
  #                              prob = 1 - out_fore_leroux$Pi[na.ind])
  
  #  # Plot the posterior means jointly across day, state:
  #  library(scales)
  #  states = 1:49
  #  dev.new(); 
  #  par(mai = c(.9,0.9,.4,0))
  #  filled.contour(x = tau, 
  #                 y = states, 
  #                 #z = t(Z_imputed),  #colMeans(out_fore$Zhat), 
  #                 z = t(sim_data$Z_full),
  #                 #z = t(Z),
  #                 #color = function(n) hcl.colors(n, "terrain"),
  #                 #color.palette=viridis_pal(), 
  #                 color.palette = terrain.colors, 
  #                 ylab = 'State ID', xlab =' Observation Point',plot.axes = { axis(2, states[seq(1, T, by = 4)], cex.axis = 1.5)
  #                   axis(1, at = round(tau[seq(1, m, by = 24)],1),  cex.axis = 1.5, mgp=c(3,0.5,0))},
  #                 cex.lab=1.5, cex.axis = 1, cex.main =1.5, 
  #                 zlim = range(Z, na.rm=TRUE),
  #                 main = 'Simulation Counts: Forecast')
  
  
  
  
  # Other estimators:
  ############# Poisson-LCAR 
  out_fore_pois = fosr.count_leroux_rho1(Z , tau = tau, X =X, W=W,
                                  K = K, use_approx_poisson = TRUE,
                                  Offset = NULL,
                                  rho.possible = seq(0, 0.95, 0.05),
                                  nsave = 5000, nburn = 5000, nskip = 4,
                                  mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                  computeDIC = TRUE)
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_pois = array(NA,c(length(id_pred),m))
  post_Zpred_pois = NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    # Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    # Z_obs =c(Z_obs,sim_data$Z_full[s,ind_fore])
    
    if(length(ind_fore)>1) 
    {Zhat_pois[d,ind_fore] = colMeans(out_fore_pois$Zhat[,s,ind_fore])}
    else {Zhat_pois[d,ind_fore] = mean(out_fore_pois$Zhat[,s,ind_fore]) }
    
    post_Zpred_pois=cbind(post_Zpred_pois, out_fore_pois$Zpred[,s, ind_fore])
  }  
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_pois$fk, 
                                              post_alpha_j = out_fore_pois$alpha[,j,])
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[2,i] =alpha_tilde_width[2,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[2,i] =alpha_tilde_prob[2,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[2,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[2,i] = alpha_tilde_width[2,i]/((p-1)*m)
  alpha_tilde_prob[2,i] =   alpha_tilde_prob[2,i]/((p-1)*m)
  
  
  
  
  ####################### Gauss-LCAR
  Y=sqrt(Z/rep(1,n))  ##
  # Y=sqrt(Z/(pop[-c(2,12),2]))  #
  
  out_fore_gauss = fosr.guass_icar(Y = Y, E=rep(1,n), tau= tau, X = X, K = K, W = W,
                                   nsave = 5000, nburn = 5000, nskip = 4,
                                   mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                   computeDIC = TRUE)
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_gauss = array(NA,c(length(id_pred),m))
  post_Zpred_gauss = NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    # Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    # Z_obs =c(Z_obs,sim_data$Z_full[s,ind_fore])
    
    if(length(ind_fore)>1) 
    {Zhat_gauss[d,ind_fore] = colMeans(out_fore_gauss$Zhat[,s,ind_fore])}
    else {Zhat_gauss[d,ind_fore] = mean(out_fore_gauss$Zhat[,s,ind_fore]) }
    
    post_Zpred_gauss=cbind(post_Zpred_gauss, out_fore_gauss$Zpred[,s, ind_fore])
  }  
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_gauss$fk, 
                                              post_alpha_j = out_fore_gauss$alpha[,j,])
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[3,i] =alpha_tilde_width[3,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[3,i] =alpha_tilde_prob[3,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[3,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[3,i] = alpha_tilde_width[3,i]/((p-1)*m)
  alpha_tilde_prob[3,i] =   alpha_tilde_prob[3,i]/((p-1)*m)
  
  
  
####################### Nb-Ind
  out_fore_ind = fosr.count_ind(Z = Z, tau = tau, X =X, K=K,
                                use_approx_poisson = FALSE, sample_r = TRUE, sample_sigma = TRUE,
                                Offset =NULL,
                                nsave = 5000, nburn = 5000, nskip = 4,
                                mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                computeDIC = TRUE)
  
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_ind = array(NA,c(length(id_pred),m))
  post_Zpred_ind = NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    # Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    # Z_obs =c(Z_obs,sim_data$Z_full[s,ind_fore])
    
    if(length(ind_fore)>1) 
    {Zhat_ind[d,ind_fore] = colMeans(out_fore_ind$Zhat[,s,ind_fore])}
    else {Zhat_ind[d,ind_fore] = mean(out_fore_ind$Zhat[,s,ind_fore]) }
    
    post_Zpred_ind=cbind(post_Zpred_ind, out_fore_ind$Zpred[,s, ind_fore])
  }  
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_ind$fk, 
                                              post_alpha_j = out_fore_ind$alpha[,j,])
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[4,i] =alpha_tilde_width[4,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[4,i] =alpha_tilde_prob[4,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[4,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[4,i] = alpha_tilde_width[4,i]/((p-1)*m)
  alpha_tilde_prob[4,i] =   alpha_tilde_prob[4,i]/((p-1)*m)
  
  
  
  

  # MCMC sampler:  
  out_fore_icar = fosr.count_icar(Z, tau, X =X, W=W,
                             K = K, use_approx_poisson = FALSE,
                             Offset = NULL,
                             nsave = 5000, nburn = 5000, nskip = 4,
                             mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "Pi", "sigma_e")
  )
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_icar = array(NA,c(length(id_pred),m))
  post_Zpred_icar = NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    #Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    #Z_obs =c(Z_obs,sim_data$Z_full[s,ind_fore])
    
    if(length(ind_fore)>1) 
    {Zhat_icar[d,ind_fore] = colMeans(out_fore_icar$Zhat[,s,ind_fore])}
    else {Zhat_icar[d,ind_fore] = mean(out_fore_icar$Zhat[,s,ind_fore]) }
    
    post_Zpred_icar = cbind(post_Zpred_icar, out_fore_icar$Zpred[,s, ind_fore])
    
  }
  
  
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_icar$fk, 
                                              post_alpha_j = out_fore_icar$alpha[,j,])
    
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
 
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[5,i] =alpha_tilde_width[5,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[5,i] =alpha_tilde_prob[5,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  
  rmse[5,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[5,i] = alpha_tilde_width[5,i]/((p-1)*m)
  alpha_tilde_prob[5,i] =   alpha_tilde_prob[5,i]/((p-1)*m)
  
  
  
  #################
  
  # MCMC sampler:  Nb-CCAR
  out_fore_cressie_rho1 = fosr.count_cressie_rho1(Z, tau, X =X, W=W,
                                                  K = K, use_approx_poisson = FALSE,
                                                  Offset = NULL,
                                                  rho.possible = seq(0, 0.95, 0.05),
                                                  nsave = 5000, nburn = 5000, nskip = 4,
                                                  mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "Pi", "sigma_e","rho")
  )
  

  # Posterior mean and distribution of the forecast:
  Zhat_cressie_rho1 = array(NA,c(length(id_pred),m))
  post_Zpred_cressie_rho1=NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    #Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    #Z_obs[d,ind_fore] =sim_data$Z_full[s,ind_fore]
    
    if(length(ind_fore)>1) 
    {Zhat_cressie_rho1[d,ind_fore] = colMeans(out_fore_cressie_rho1$Zhat[,s,ind_fore])}
    else {Zhat_cressie_rho1[d,ind_fore] = mean(out_fore_cressie_rho1$Zhat[,s,ind_fore]) }
    
    post_Zpred_cressie_rho1= cbind(post_Zpred_cressie_rho1,out_fore_cressie_rho1$Zpred[,s, ind_fore])
  }
 

  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_cressie_rho1$fk, 
                                              post_alpha_j = out_fore_cressie_rho1$alpha[,j,])
    
    # Effective size:
    print(getEffSize(post_alpha_tilde_j))
    
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
 
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[6,i] =alpha_tilde_width[6,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[6,i] =alpha_tilde_prob[6,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[6,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[6,i] = alpha_tilde_width[6,i]/((p-1)*m)
  alpha_tilde_prob[6,i] =   alpha_tilde_prob[6,i]/((p-1)*m)
  
  
  
  
  library(refund)
  ############# FPCA-basis
  Y=sqrt(Z/rep(1, n))  ##
  out_fore_fpc = fosr.guass_fpc(Y=Y, E=rep(1, n), tau= tau, X = X, K = K, W= W, pve=0.99,
                                nsave = 5000, nburn = 5000, nskip = 4,
                                mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                computeDIC = TRUE)
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_fpc = array(NA,c(length(id_pred),m))
  post_Zpred_fpc=NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    #Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    #Z_obs[d,ind_fore] =sim_data$Z_full[s,ind_fore]
    
    if(length(ind_fore)>1) 
    {Zhat_fpc[d,ind_fore] = colMeans(out_fore_fpc$Zhat[,s,ind_fore])}
    else {Zhat_fpc[d,ind_fore] = mean(out_fore_fpc$Zhat[,s,ind_fore]) }
    
    post_Zpred_fpc= cbind(post_Zpred_fpc,out_fore_fpc$Zpred[,s, ind_fore])
  }
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_fpc$fk, 
                                              post_alpha_j = out_fore_fpc$alpha[,j,])
    
    
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[7,i] =alpha_tilde_width[7,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[7,i] =alpha_tilde_prob[7,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[7,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[7,i] = alpha_tilde_width[7,i]/((p-1)*m)
  alpha_tilde_prob[7,i] =   alpha_tilde_prob[7,i]/((p-1)*m)
  
  
  
  
  ############## spline-basis
  out_fore_spline = fosr.count_icar_spline(Z = Z, tau = tau, X =X, W=W,
                                           use_approx_poisson = FALSE, sample_r = TRUE, sample_sigma = TRUE,
                                           Offset = NULL,
                                           nsave = 5000, nburn = 5000, nskip = 4,
                                           mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                           computeDIC = TRUE)
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_spline = array(NA,c(length(id_pred),m))
  post_Zpred_spline=NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    #Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    #Z_obs[d,ind_fore] =sim_data$Z_full[s,ind_fore]
    
    if(length(ind_fore)>1) 
    {Zhat_spline[d,ind_fore] = colMeans(out_fore_spline$Zhat[,s,ind_fore])}
    else {Zhat_spline[d,ind_fore] = mean(out_fore_spline$Zhat[,s,ind_fore]) }
    
    post_Zpred_spline= cbind(post_Zpred_spline,out_fore_spline$Zpred[,s, ind_fore])
  }
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_spline$fk, 
                                              post_alpha_j = out_fore_spline$alpha[,j,])
    
    
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[8,i] =alpha_tilde_width[8,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[8,i] =alpha_tilde_prob[8,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[8,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[8,i] = alpha_tilde_width[8,i]/((p-1)*m)
  alpha_tilde_prob[8,i] =   alpha_tilde_prob[8,i]/((p-1)*m)
  
  
  
  ############## spline-NIG prior
  out_fore_spline_nig = fosr.count_icar_spline_nig(Z = Z, tau = tau, X =X, W=W,
                                                   use_approx_poisson = FALSE, sample_r = TRUE, sample_sigma = TRUE,
                                                   Offset =NULL,
                                                   nsave = 5000, nburn = 5000, nskip = 4,
                                                   mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                                   computeDIC = TRUE)
  
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_spline_nig = array(NA,c(length(id_pred),m))
  post_Zpred_spline_nig=NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    #Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    #Z_obs[d,ind_fore] =sim_data$Z_full[s,ind_fore]
    
    if(length(ind_fore)>1) 
    {Zhat_spline_nig[d,ind_fore] = colMeans(out_fore_spline_nig$Zhat[,s,ind_fore])}
    else {Zhat_spline_nig[d,ind_fore] = mean(out_fore_spline_nig$Zhat[,s,ind_fore]) }
    
    post_Zpred_spline_nig= cbind(post_Zpred_spline_nig,out_fore_spline_nig$Zpred[,s, ind_fore])
  }
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_spline_nig$fk, 
                                              post_alpha_j = out_fore_spline_nig$alpha[,j,])
    
    
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[9,i] =alpha_tilde_width[9,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[9,i] =alpha_tilde_prob[9,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[9,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[9,i] = alpha_tilde_width[9,i]/((p-1)*m)
  alpha_tilde_prob[9,i] =   alpha_tilde_prob[9,i]/((p-1)*m)
  
  
  
  ##################  NIG
  out_fore_icar_nig = fosr.count_icar_nig(Z = Z, tau = tau, X =X, K=K, W=W,
                                          use_approx_poisson = FALSE, sample_r = TRUE, sample_sigma = TRUE,
                                          Offset = NULL,
                                          nsave = 5000, nburn = 5000, nskip = 4,
                                          mcmc_params = list("beta","fk", "alpha", "Zhat", "Zpred", "r", "sigma_e"),
                                          computeDIC = TRUE)
  
  
  
  
  # Posterior mean and distribution of the forecast:
  Zhat_icar_nig = array(NA,c(length(id_pred),m))
  post_Zpred_icar_nig=NULL
  
  for(d in 1:length(id_pred)){
    s=id_pred[d]
    ind_fore=(m_fore-date[s]+2):m    # indices to forecast
    #Z_true[d,ind_fore]=sim_data$Z_true[s,ind_fore]
    #Z_obs[d,ind_fore] =sim_data$Z_full[s,ind_fore]
    
    if(length(ind_fore)>1) 
    {Zhat_icar_nig[d,ind_fore] = colMeans(out_fore_icar_nig$Zhat[,s,ind_fore])}
    else {Zhat_icar_nig[d,ind_fore] = mean(out_fore_icar_nig$Zhat[,s,ind_fore]) }
    
    post_Zpred_icar_nig= cbind(post_Zpred_icar_nig,out_fore_icar_nig$Zpred[,s, ind_fore])
  }
  
  ############ alpha_j(t)
  # Plot the regression coefficient functions:
  alpha_tilde_hat = matrix(0, m,p-1)  ## store regression coefficient functions 
  
  for(j in 2:ncol(X)){ # skip the intercept
    # Obtain the posterior draws for the (non-dynamic) function, predictor j:
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out_fore_icar_nig$fk, 
                                              post_alpha_j = out_fore_icar_nig$alpha[,j,])
    
    
    # Plot the posterior mean and credible intervals/bands:
    # Pointwise intervals:
    ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
    # Pointwise bands:
    cb_j = credBands(post_alpha_tilde_j)
    
    
    alpha_tilde_hat[,j-1] = colMeans(post_alpha_tilde_j)
    
    alpha_tilde_width[10,i] =alpha_tilde_width[10,i] + sum(abs(ci_j[,2] - ci_j[,1]))
    
    alpha_tilde_prob[10,i] =alpha_tilde_prob[10,i] + sum((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
  }
  
  rmse[10,i] = sqrt(mean((alpha_tilde_hat - sim_data$alpha_tilde_true[,-1])^2))
  alpha_tilde_width[10,i] = alpha_tilde_width[10,i]/((p-1)*m)
  alpha_tilde_prob[10,i] =   alpha_tilde_prob[10,i]/((p-1)*m)
  
  
  
  
  
  ###############
  # Store the errors:
  errors_baye[1,i] = mean(abs(Z_true - Zhat_leroux_rho1),na.rm=TRUE)
  errors_baye[2,i] = mean(abs(Z_true - Zhat_pois),na.rm=TRUE)
  errors_baye[3,i] = mean(abs(Z_true - Zhat_gauss),na.rm=TRUE)
  errors_baye[4,i] = mean(abs(Z_true - Zhat_ind),na.rm=TRUE)
 
 
  errors_Y[1,i] = mean(abs(Z_obs - apply(post_Zpred_leroux_rho1,2,median)))   ## 后验中位数作为Yhat，再计算MAE
  errors_Y[2,i] = mean(abs(Z_obs - apply(post_Zpred_pois,2,median))) 
  errors_Y[3,i] = mean(abs(Z_obs - apply(post_Zpred_gauss,2,median))) 
  errors_Y[4,i] = mean(abs(Z_obs - apply(post_Zpred_ind,2,median))) 
  
  
  # Intervals:
  ci_nb = HPDinterval(as.mcmc(post_Zpred_leroux_rho1))
  ci_pois = HPDinterval(as.mcmc(post_Zpred_pois))
  ci_gauss = HPDinterval(as.mcmc(post_Zpred_gauss))
  ci_ind = HPDinterval(as.mcmc(post_Zpred_ind))

  
  # Coverage prob:
  cover_prob_baye[1, i] =   mean((Z_obs >= ci_nb[,1])*(Z_obs <= ci_nb[,2]))
  cover_prob_baye[2, i] =   mean((Z_obs >= ci_pois[,1])*(Z_obs <= ci_pois[,2]))
  cover_prob_baye[3, i] =   mean((Z_obs >= ci_gauss[,1])*(Z_obs <= ci_gauss[,2]))
  cover_prob_baye[4, i] =   mean((Z_obs >= ci_ind[,1])*(Z_obs <= ci_ind[,2]))

  
  # Coverage widths: 
  cover_width_baye[1, i] = median(abs(ci_nb[,2] - ci_nb[,1]))
  cover_width_baye[2, i] = median(abs(ci_pois[,2] - ci_pois[,1]))
  cover_width_baye[3, i] = median(abs(ci_gauss[,2] - ci_gauss[,1]))
  cover_width_baye[4, i] = median(abs(ci_ind[,2] - ci_ind[,1]))


  
  
}



