rm(list = ls())
#------------------------------------------------------------------------------------------------------------------
# Model specifications:
#use_approx_poisson = FALSE # Poisson or Negative-Binomial?
K = 8 # Number of (unknown) basis functions 
#------------------------------------------------------------------------------------------------------------------
# Load the libraries:
library(dfosr); # devtools::install_github("drkowal/dfosr")
library(coda) # For MCMC output
library(quantmod) # For population data (from FRED)
library(fds); library(viridis) # For plotting
library(ggplot2)
library(forecast)
library(tscount)
library(refund)

library(KFAS)
library(truncdist)
library(BayesLogit)

#setwd("D://R code_integer valued paper/dfosr.nb")
#devtools::load_all(".") #library(dfosr.nb)

# Read the data:

daily_cases2=readxl::read_excel("D:/daily-cases.xlsx",2)

######covariates
X0 = readxl::read_excel("D:/covariates.xlsx")
colnames(X0)=c("states","Population Density","Percentage of People over 65","State-level GDP","Transit change")
X49=cbind(1,scale(X0[,2:5]))  ## standardized

# Population:
pop = read.csv("D:/states_population.csv")


m=350
m_fore=375

###complete cases
cases=as.matrix(daily_cases2[-c(2,12),2:(m+1)])  ### 分析天数m
####
Z=cases

tau = c(1:m)  ## time period
T = n = 49    ##  
states = 1:49

############################## 
date_covid = readxl::read_excel("D:/case_date.xlsx")
date = as.numeric(as.matrix(date_covid[,2]))


####   Now delete the observations:
id_pred=NULL  ### id to be predicted
for(ii in 1:n){
  if(date[ii]>(m_fore-m+1))
  { id_pred=rbind(id_pred,ii)
  Z[ii, (m_fore-date[ii]+2):m ]=NA
  }
}



###### 
par(mar=c(4.5,4.5,3,1))    
plot(cases[1,]/100,ylim=c(0,600),lwd=1,type='l',col="gray", xlab='Days since the first 20 cases', cex.lab=1.3, cex.main=1.5,
     ylab=expression(paste('Daily Increased Cases (hundreds)')),main='COVID-19 Case Counts in the States')
for(i in 1:49){
  if(i != c(1,44,33,31,10,3,5,2,12) )
    lines(cases[i,]/100,lwd=1,type='l',col="gray")
}
lines(cases[42,]/100,lwd=2,type='l',col='plum2')
lines(cases[31,]/100,lwd=2,type='l',col=4)
lines(cases[29,]/100,lwd=2,type='l',col=rainbow(49)[25])
lines(cases[9,]/100,lwd=2,type='l',col=3)
lines(cases[2,]/100,lwd=2,type='l',col='chocolate1')
lines(cases[4,]/100,lwd=2,type='l',col='darkgoldenrod1' )

legend('topleft',ncol=2,legend = c("Texas",'New York','New Jersey','Florida','Arizona','California','Other states'), col=c('plum2',4, rainbow(49)[25],3, 'chocolate1','darkgoldenrod1','gray'),
       lwd=2, x.intersp=0.4, text.width=35,seg.len=1)



################## 邻接矩阵W
library(maps)
library(maptools)
library(spdep)
library(classInt) ## Will be used for plotting maps later
library(RColorBrewer) ## Will be used for colored maps

usa.state = map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"), function(x) x[1])   ##63
usa.poly = map2SpatialPolygons(usa.state, IDs=state.ID)    ##49
usa.nb = poly2nb(usa.poly)
W = usa.adj.mat = nb2mat(usa.nb, style="B") ### 邻接矩阵（49个州，不含AK,HI）
write.csv(usa.adj.mat,"D:/49states_adjmatrix.csv")


W=read.csv("D:/49states_adjmatrix.csv")
statesname=W[,1]
W = as.matrix(W[,-1])

########## 预装函数
setwd("D:/fosr.count_spatial")
source("helper_functions_spatial.R")
source("mcmc_sampler_spatial.R")


#------------------------------------------------------------------------------------------------------------------
# Run the MCMC: count versions
#------------------------------------------------------------------------------------------------------------------
#set.seed(1234) # For reproducibility

###
out = fosr.count_leroux_rho1(Z = Z, tau = tau, X =X49, W=W,
                              K = K, use_approx_poisson = FALSE,
                              Offset = matrix(rep(as.matrix(pop[-c(2,12),2]), times = m), nr = T),
                              rho.possible = seq(0, 0.95, 0.05),
                              nsave = 50, nburn = 10, nskip = 4,
                              mcmc_params = list("beta","fk", "alpha","mu_k", "Zhat", "Zpred", "r", "sigma_e","rho"), 
                              computeDIC = TRUE)


# Posterior mean and distribution of the forecast:
Zhat = array(NA,c(length(id_pred),m))
post_Zpred =NULL

for(d in 1:length(id_pred)){
  s=id_pred[d]
  ind_fore=(m_fore-date[s]+2):m    # indices to forecast

  
  if(length(ind_fore)>1) 
  {Zhat[d,ind_fore] = colMeans(out$Zhat[,s,ind_fore])}   ## point estimation
  else {Zhat[d,ind_fore] = mean(out$Zhat[,s,ind_fore]) }
  
  post_Zpred = cbind(post_Zpred, out$Zpred[,s, ind_fore])   ## posterior samples, for interval estimation

}


s=14
ind_fore=(m_fore-date[s]+2):m  

# Summary plot:
plot_curve(out$Zpred[,s,], tau, include_joint = FALSE, main = 'Iowa'); 
lines(tau, Z[s,], type='p', pch = 1); 
abline(v = tau[ind_fore[1]-1], lty = 6); 
lines(tau[ind_fore], cases[s, ind_fore], type ='p',pch=2)



############### 

# Plot the regression coefficient functions:
for(j in 2:ncol(X49)){ # skip the intercept
  # Obtain the posterior draws for the (non-dynamic) function, predictor j:
  post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out$fk, 
                                            post_alpha_j = out$alpha[,j,])
  
  
  # Plot the posterior mean and credible intervals/bands:
  # Pointwise intervals:
  ci_j = t(apply(post_alpha_tilde_j, 2, quantile, c(0.05/2, 1 - 0.05/2)))
  # Pointwise bands:
  cb_j = credBands(post_alpha_tilde_j)
  
  dev.new(); par(mai = c(.8,.9,.4,.4));
  plot(tau, post_alpha_tilde_j[1,], type = "n", ylim = range(ci_j, cb_j, 0, na.rm = TRUE),
       xlab = 't', ylab = "", main = colnames(X49)[j],
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  #axis(1, at = mon.ind, labels = substr(month.name, 0, 3), cex.axis = 2)
  polygon(c(tau, rev(tau)), c(cb_j[, 2], rev(cb_j[, 1])), col = "gray50",border = NA)
  polygon(c(tau, rev(tau)), c(ci_j[, 2], rev(ci_j[, 1])), col = "gray",border = NA)
  lines(tau, colMeans(post_alpha_tilde_j), lwd = 3, col = 'chocolate1')
  abline(h = 0, lty = 3, lwd=3)
  
}

# Plot the factors:
plot_factors(post_beta = out$beta, states)   

# Plot the loading curves:
plot_flc(post_fk = out$fk, tau = tau)
