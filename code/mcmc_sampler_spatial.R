
##########################
fosr.count_leroux_rho1 = function(Z, tau, X = NULL, K = NULL, W,
                                  use_approx_poisson = TRUE, sample_r = TRUE, sample_sigma = TRUE,
                                  Offset = NULL,
                                  rho.possible = rho.possible,
                                  nsave = 1000, nburn = 1000, nskip = 2,
                                  mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred","rho"),
                                  computeDIC = TRUE){
  
  # Some options (for now):
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?
  
  
  ############# rho
  n.rho <- length(rho.possible)
  mat.Q <- rep(0, length(rho.possible))
  det.Q <- rep(0, length(rho.possible))
  mat.Q <- as.list(mat.Q)
  
  W.rowsum <- apply(W, 1, sum)
  
  #### Calculate the possible values for the precision matrix Q and |Q|
  #####################################################################
  
  for(h in 1:length(rho.possible))
  {
    Q <- -rho.possible[h] * W
    diag(Q) <- W.rowsum * rho.possible[h] + 1 - rho.possible[h]
    mat.Q[[h]] <- Q
    eigen.temp <- eigen(Q)
    det.Q[h] <- 0.5 * sum(log(eigen.temp$values))
  } 
  
  
  
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  
  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    #Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100
  
  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  
  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value
  
  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi
  
  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)
  
  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1
  
  # Impute the missing points:
  if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                      size = r,
                                      prob = 1 - Pi[na.ind])
  
  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  
  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals
  
  rho=rep(0.6, 1)
  
  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef
    
    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }
  
  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])
  
  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)
  
  
  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K);
  sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
  
  # Update the error SD for eta / conditional SD for gamma:
  sigma_gamma_tk = matrix(rep(sigma_delta_k, each = T)/
                            sqrt( rep(rho*rowSums(W)+1-rho,K)), nrow = T)
  
  #---------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
    
    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term
    
    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term
    
    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('rho', mcmc_params)) || computeDIC) post.rho = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))
    
    # Transform for FDLM parameters:
    BtTheta = tcrossprod(t(splineInfo$Bmat), Theta - logOffset)
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtTheta,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_e^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;
    
    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtTheta, Psi);    sigma_tilde = sigma_e
    
    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 )   ###修改
      gamma_k = gamma_tk[,k]
      
      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        #Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }
    
    # And sample the errors gamma_tk:
    
    ######
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_tk^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2
    
    for(t in 1:T){
      ell_gamma_t = ell_gamma_tk1[t,] + rho*(W %*% gamma_tk)[t,]/sigma_delta_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_t * postSD[t,]^2   ##逐个元素相乘 只有一行
      
      gamma_tk[t, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[t,]))  ##按行更新，迭代新值进入下一行的更新
    }
    
    
    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk
    
    ##update eta
    eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(rho/(rowSums(W)* rho +1-rho),K), T)
    
    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)
      
    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]
    
    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }
    
    # Variance part:   
    # Standardize, then reconstruct as matrix of size T x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(eta_tk*sqrt(rep((rho*rowSums(W)+1-rho),K)), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }
    
    # Update the error SD for gamma:
    sigma_gamma_tk = matrix(rep(sigma_delta_k, each = T)/
                              sqrt( rep(rho*rowSums(W)+1-rho,K)), nrow = T)
    
    #### Sample from rho
    ####################
    
    log.rho.prob <- rep(0, n.rho)
    for(j in 1:n.rho)
    { 
      log.rho.prob[j] = K* det.Q[j]
      for(k in 1:K)
      {
        log.rho.prob[j] <- log.rho.prob[j] - 0.5 * t(gamma_tk[,k]) %*% mat.Q[[j]] %*% gamma_tk[,k] / (sigma_delta_k[k])^2
        
      }
    }
    
    rho.prob <- exp(log.rho.prob) / sum(exp(log.rho.prob))
    rho <- sample(x=rho.possible, size=1, prob=rho.prob)
    
    
    
    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      
      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
      
      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('rho', mcmc_params)) || computeDIC) post.rho[isave,] = rho
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred
  if(!is.na(match('rho', mcmc_params))) mcmc_output$rho = post.rho
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}




###################

fosr.count_cressie_rho1 = function(Z, tau, X = NULL, K = NULL, W,
                                   use_approx_poisson = TRUE, sample_r = TRUE, sample_sigma = TRUE,
                                   Offset = NULL,
                                   rho.possible = rho.possible,
                                   nsave = 1000, nburn = 1000, nskip = 2,
                                   mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred","rho"),
                                   computeDIC = TRUE){
  
  # Some options (for now):
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?
  
  ############# rho
  n.rho <- length(rho.possible)
  mat.Q <- rep(0, length(rho.possible))
  det.Q <- rep(0, length(rho.possible))
  mat.Q <- as.list(mat.Q)
  
  W.rowsum <- apply(W, 1, sum)
  
  #### Calculate the possible values for the precision matrix Q and |Q|
  #####################################################################
  
  for(h in 1:length(rho.possible))
  {
    Q <- -rho.possible[h] * W
    diag(Q) <- W.rowsum
    mat.Q[[h]] <- Q
    eigen.temp <- eigen(Q)
    det.Q[h] <- 0.5 * sum(log(eigen.temp$values))
  } 
  
  rho=rep(0.8, 1)
  
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  
  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    #Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100
  
  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  
  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value
  
  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi
  
  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)
  
  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1
  
  # Impute the missing points:
  if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                      size = r,
                                      prob = 1 - Pi[na.ind])
  
  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  
  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals
  
  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef
    
    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }
  
  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])
  
  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)
  
  
  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K);
  sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
  
  # Update the error SD for eta / conditional SD for gamma:
  sigma_gamma_tk = matrix(rep(sigma_delta_k, each = T)/sqrt(rep(rowSums(W), K)), nrow = T)
  
  #---------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
    
    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term
    
    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term
    
    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('rho', mcmc_params)) || computeDIC) post.rho = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))
    
    # Transform for FDLM parameters:
    BtTheta = tcrossprod(t(splineInfo$Bmat), Theta - logOffset)
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtTheta,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_e^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;
    
    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtTheta, Psi); sigma_tilde = sigma_e
    
    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 )   ###修改
      gamma_k = gamma_tk[,k]
      
      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }
    
    # And sample the errors gamma_tk:
    
    ######
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_tk^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2
    
    for(t in 1:T){
      ell_gamma_t = ell_gamma_tk1[t,] + rho*(W %*% gamma_tk)[t,]/sigma_delta_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_t * postSD[t,]^2   ##逐个元素相乘 只有一行
      
      gamma_tk[t, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[t,]))  ##按行更新，迭代新值进入下一行的更新
    }
    
    
    
    
    
    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk
    
    ##update eta
    eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(rho/rowSums(W),K), T)
    
    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)
      
    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]
    
    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }
    
    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(eta_tk*sqrt(rep(rowSums(W),K)), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }
    
    
    
    # Update the error SD for gamma:
    sigma_gamma_tk = rep(sigma_delta_k, each = T)/sqrt( matrix(rep(rowSums(W),K),nrow=T) )
    
    
    #### Sample from rho
    ####################
    
    log.rho.prob <- rep(0, n.rho)
    for(j in 1:n.rho)
    { 
      log.rho.prob[j] = K* det.Q[j]
      for(k in 1:K)
      {
        log.rho.prob[j] <- log.rho.prob[j] - 0.5 * t(gamma_tk[,k]) %*% mat.Q[[j]] %*% gamma_tk[,k] / (sigma_delta_k[k])^2
        
      }
    }
    
    rho.prob <- exp(log.rho.prob) / sum(exp(log.rho.prob))
    rho <- sample(x=rho.possible, size=1, prob=rho.prob)
    
    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      
      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
      
      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('rho', mcmc_params)) || computeDIC) post.rho[isave,] = rho
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred
  if(!is.na(match('rho', mcmc_params))) mcmc_output$rho = post.rho
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}










#########################
#' MCMC Sampling Algorithm for the Count-valued Function-on-Scalars Regression Model
#'
#' Runs the MCMC for the function-on-scalars regression model with count responses based on
#' a reduced-rank expansion. Here we assume the factor regression has ICAR errors,
#' as well as some additional default conditions.
#'
#' @param Z the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param use_approx_poisson logical; if TRUE, set the size parameter \code{r = 1000} to approximate the Poisson distribution
#' @param sample_r logical; if TRUE, sample the size parameter \code{r},
#' otherwise compute a method-of-moments estimator (unless \code{use_approx_poisson = TRUE})
#' @param sample_sigma logical; if TRUE, sample the variance parameter \code{sigma},
#' otherwise use \code{sigma = 0.1} to include a small amount of small-scale variability
#' @param Offset the \code{T x m} matrix of offsets; if NULL, assume no offsets
#' @param W adjacent matrix
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept)
#' \item "r" (Negative Binomial size parameter)
#' \item "sigma_e" (factor-level error standard deviation)
#' \item "Zhat" (conditional expectation of the count observations, Z)
#' \item "Zpred" (posterior predictive distribution of Z)
#' \item "Zfore" (one-step forecast)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @examples
#' \dontrun{
#' # Fixme
#'}
#'
#' @import KFAS truncdist dfosr
#' @importFrom BayesLogit rpg
# @export
fosr.count_icar = function(Z, tau, X = NULL, K = NULL, W,
                           use_approx_poisson = TRUE, sample_r = TRUE, sample_sigma = TRUE,
                           Offset = NULL,
                           nsave = 1000, nburn = 1000, nskip = 2,
                           mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                           computeDIC = TRUE){

  # Some options (for now):
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?

  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    #Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100

  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)

  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1

  # Impute the missing points:
  if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                      size = r,
                                      prob = 1 - Pi[na.ind])

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)

  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)


  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K);
  sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))

  # Update the error SD for eta / conditional SD for gamma:
  sigma_gamma_tk = matrix(rep(sigma_delta_k, each = T)/sqrt(rep(rowSums(W), K)), nrow = T)

  #---------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term

    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term

    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))

    # Transform for FDLM parameters:
    BtTheta = tcrossprod(t(splineInfo$Bmat), Theta - logOffset)
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtTheta,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_e^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtTheta, Psi); sigma_tilde = sigma_e

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 )   ###修改
      gamma_k = gamma_tk[,k]

      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_tk:

    # postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_tk^2))    ###vector
    # postMean = matrix((Y_tilde - X%*%alpha_pk)/rep(sigma_tilde^2, times = K))*postSD^2   ##vetor form
    # gamma_tk = matrix(rnorm(n = T*K, mean = postMean, sd = postSD), nrow = T)


    ######
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_tk^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2

    for(t in 1:T){
      ell_gamma_t = ell_gamma_tk1[t,] + (W %*% gamma_tk)[t,]/sigma_delta_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_t * postSD[t,]^2   ##逐个元素相乘 只有一行

      gamma_tk[t, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[t,]))  ##按行更新，迭代新值进入下一行的更新
    }


    
    
    
    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk

    ##update eta
    eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)

    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)

    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(eta_tk*sqrt(rep(rowSums(W),K)), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }



    # Update the error SD for gamma:
    sigma_gamma_tk = rep(sigma_delta_k, each = T)/sqrt( matrix(rep(rowSums(W),K),nrow=T) )
    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}







################ Gauss-ICAR 整数变换为实数

#' MCMC Sampling Algorithm for the Function-on-Scalars Regression Model
#'
#' Runs the MCMC for the function-on-scalars regression model based on
#' an FDLM-type expansion. Here we assume the factor regression has ICAR errors,
#' which allows for subject-specific random effects,
#' as well as some additional default conditions.
#'
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of subjects and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param E the vector of the population in states
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{n x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param W adjacent matrix
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept term for factor k)
#' \item "sigma_e" (observation error SD)
#' \item "sigma_g" (random effects SD)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive distribution of Y)
#' \item "Zhat" (fitted values)
#' \item "Zpred" (posterior predictive distribution of Z)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{nm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x n x M},  may be inefficient
#'
#' @examples
#' # Simulate some data:
#' sim_data = simulate_fosr(n = 100, m = 20, p_0 = 100, p_1 = 5)
#'
#' # Data:
#' Y = sim_data$Y; X = sim_data$X; tau = sim_data$tau
#'
#' # Dimensions:
#' n = nrow(Y); m = ncol(Y); p = ncol(X)
#'
#' # Run the FOSR:
#' out = fosr(Y = Y, tau = tau, X = X, K = 6, mcmc_params = list("fk", "alpha", "Yhat"))
#'
#' # Plot a posterior summary of a regression function, say j = 3:
#' j = 3; post_alpha_tilde_j = get_post_alpha_tilde(out$fk, out$alpha[,j,])
#' plot_curve(post_alpha_tilde_j, tau = tau)
#' # Add the true curve:
#' lines(tau, sim_data$alpha_tilde_true[,j], lwd=6, col='green', lty=6)
#'
#' # Plot the loading curves:
#' plot_flc(out$fk, tau = tau)
#'
#' # Plot the fitted values for a random subject:
#' i = sample(1:n, 1)
#' plot_fitted(y = Y[i,], mu = colMeans(out$Yhat[,i,]),
#'             postY = out$Yhat[,i,], y_true = sim_data$Y_true[i,], t01 = tau)
#'
#' @import KFAS truncdist dfosr
#' @importFrom BayesLogit rpg
#' @export
fosr.guass_icar = function(Y, E, tau, X = NULL, K = NULL, W,
                nsave = 1000, nburn = 1000, nskip = 3,
                mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                computeDIC = TRUE){

  # Some options (for now):
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?

  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  n = nrow(Y); m = ncol(Y); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms:

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing){na.ind = which(is.na(Yna), arr.ind = TRUE); Y = inits$Y0}
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  Yhat = tcrossprod(Beta, Fmat)

  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, n), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)

  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_ik = matrix(0, nrow = n, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= n) or OLS (p < n)
  for(k in 1:K) {
    if(p >= n){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_et,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Y, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_ik[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # Initialize the corresponding SD term(s):

  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_ik = gamma_ik - (W %*% gamma_ik) * matrix(rep(1/rowSums(W),K), n)


  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))

  # Update the error SD for gamma:
  sigma_gamma_ik = matrix(rep(sigma_delta_k, each = n)/sqrt(rep(rowSums(W), K)), nrow = n)

  #----------------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term

    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term

    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g = array(NA, c(nsave, n, K))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Ypred', mcmc_params))|| computeDIC) post.Ypred = array(NA, c(nsave, n, m))
  if(!is.na(match('Zhat', mcmc_params))|| computeDIC) post.Zhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Zpred', mcmc_params))|| computeDIC) post.Zpred = array(NA, c(nsave, n, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    if(any.missing){
      Y[na.ind] = Yhat[na.ind] + sigma_et[na.ind[,1]]*rnorm(nrow(na.ind))
      BtY = tcrossprod(t(splineInfo$Bmat), Y)
    }
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    #tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k];
      sigma_tilde_k = sqrt(sigma_tilde^2 )
      gamma_k = gamma_ik[,k]

      if(p >= n){
        # Fast sampler for p >= n (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
       } else {
      # Fast sampler for p < T (Rue, 2001?)
      if(p > 1){
        chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
      } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
      ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
      alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
    }
  }


    # And sample the errors gamma_ik:
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_ik^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2


    for(i in 1:n){
      ell_gamma_i = ell_gamma_tk1[i,] + (W %*% gamma_ik)[i,]/sigma_delta_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_i * postSD[i,]^2   ##逐个元素相乘 只有一行

      gamma_ik[i, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[i,]))  ##按行更新，迭代新值进入下一行的更新
    }


    # Update the factors:
    Beta = X%*%alpha_pk + gamma_ik

    ##update eta
    eta_ik = gamma_ik - (W %*% gamma_ik) * matrix(rep(1/rowSums(W),K), n)

    # And the fitted curves:
    Yhat = tcrossprod(Beta, Fmat)
    Zhat = Yhat^2 * E
    #----------------------------------------------------------------------------
    # Step 4: Sample the observation error variance
    #----------------------------------------------------------------------------
    # Or use uniform prior?
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)

    #----------------------------------------------------------------------------
    # Step 5: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Variance part:
    # Standardize, then reconstruct as matrix of size n x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(eta_ik*sqrt(rep(rowSums(W),K)), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))


    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }


    # Update the error SD for gamma:
    sigma_gamma_ik = rep(sigma_delta_k, each = n)/sqrt( matrix(rep(rowSums(W),K),nrow=n) )
    #----------------------------------------------------------------------------
    # Step 6: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 7: Adjust the ordering
    #----------------------------------------------------------------------------
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g[isave,,] = sigma_gamma_ik
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat
        if(!is.na(match('Ypred', mcmc_params)) || computeDIC) post.Ypred[isave,,] = Yhat + sigma_e*rnorm(length(Y))
        if(!is.na(match('Zhat', mcmc_params)) || computeDIC) post.Zhat[isave,,] = Zhat
        if(!is.na(match('Zpred', mcmc_params)) || computeDIC) post.Zpred[isave,,] = (Yhat + sigma_e*rnorm(length(Y)))^2 * E
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Yhat), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('sigma_g', mcmc_params))) mcmc_output$sigma_g = post.sigma_g
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_e), m*n),
                            log = TRUE), na.rm=TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}






################ Gauss-FPCA-ICAR 整数变换为实数, fpca作为基函数

#' MCMC Sampling Algorithm for the Function-on-Scalars Regression Model
#' with pre-fixed basis functions
#'
#' Runs the MCMC for the function-on-scalars regression model based on
#' a known basis expansion(functional principal components). Here we assume the factor regression has ICAR errors,
#' which allows for subject-specific random effects,
#' as well as some additional default conditions.
#'
#'
#'
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of subjects and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param E the vector of the population in states
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{n x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of FPCs to use; if NULL, select to explain \code{pve}% of variability
#' (ignored for other basis functions)
#' @param pve proportion of variability explained for the FPC basis; only used if \code{K} is NULL
#' @param W adjacent matrix
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept term for factor k)
#' \item "sigma_e" (observation error SD)
#' \item "Zhat" (fitted values)
#' \item "Zpred" (posterior predictive distribution of Z)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{nm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x n x M},  may be inefficient
#'
#' @examples
#'
#' @import KFAS truncdist dfosr
#' @importFrom refund fpca.face
#' @importFrom BayesLogit rpg
#' @export
fosr.guass_fpc = function(Y, E, tau, X = NULL, K = NULL, W, pve,
                           nsave = 1000, nburn = 1000, nskip = 3,
                           mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                           computeDIC = TRUE){

  # Some options (for now):
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?

  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  n = nrow(Y); m = ncol(Y); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms:

  ##### fpca basis
  fpca = fpca.face(Y, center = TRUE,
                   argvals = as.numeric(tau01),
                   knots=  max(20, min(ceiling(floor(median(rowSums(!is.na(Y))))/4), 150)),
                   pve=0.99)

  Fmat = fpca$efunctions
  Y0 = fpca$Yhat   # imputed Y
  # Initialize the factors

  Beta = tcrossprod(Y0, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value

  # Initialize the conditional expectation:
  Yhat = tcrossprod(Beta, Fmat)

  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing){na.ind = which(is.na(Yna), arr.ind = TRUE) }


 #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, n), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)

  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_ik = matrix(0, nrow = n, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= n) or OLS (p < n)
  for(k in 1:K) {
    if(p >= n){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_et,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Y, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_ik[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # Initialize the corresponding SD term(s):

  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_ik = gamma_ik - (W %*% gamma_ik) * matrix(rep(1/rowSums(W),K), n)


  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))

  # Update the error SD for gamma:
  sigma_gamma_ik = matrix(rep(sigma_delta_k, each = n)/sqrt(rep(rowSums(W), K)), nrow = n)

  #----------------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term

    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term

    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g = array(NA, c(nsave, n, K))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Ypred', mcmc_params)) || computeDIC ) post.Ypred = array(NA, c(nsave, n, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, n, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    if(any.missing){
      Y[na.ind] = Yhat[na.ind] + sigma_et[na.ind[,1]]*rnorm(nrow(na.ind))
      }
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:

    Y_tilde = tcrossprod(Y, t(Fmat))
    sigma_tilde = sigma_et

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k];
      sigma_tilde_k = sqrt(sigma_tilde^2 )
      gamma_k = gamma_ik[,k]

      if(p >= n){
        # Fast sampler for p >= n (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }


    # And sample the errors gamma_ik:
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_ik^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2


    for(i in 1:n){
      ell_gamma_i = ell_gamma_tk1[i,] + (W %*% gamma_ik)[i,]/sigma_delta_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_i * postSD[i,]^2   ##逐个元素相乘 只有一行

      gamma_ik[i, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[i,]))  ##按行更新，迭代新值进入下一行的更新
    }


    # Update the factors:
    Beta = X%*%alpha_pk + gamma_ik

    ##update eta
    eta_ik = gamma_ik - (W %*% gamma_ik) * matrix(rep(1/rowSums(W),K), n)

    # And the fitted curves:
    Yhat = tcrossprod(Beta, Fmat)
    Zhat = Yhat^2 * E
    #----------------------------------------------------------------------------
    # Step 4: Sample the observation error variance
    #----------------------------------------------------------------------------
    # Or use uniform prior?
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)

    #----------------------------------------------------------------------------
    # Step 5: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Variance part:
    # Standardize, then reconstruct as matrix of size n x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(eta_ik*sqrt(rep(rowSums(W),K)), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))


    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }


    # Update the error SD for gamma:
    sigma_gamma_ik = rep(sigma_delta_k, each = n)/sqrt( matrix(rep(rowSums(W),K),nrow=n) )
    #----------------------------------------------------------------------------
    # Step 6: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 7: Adjust the ordering
    #----------------------------------------------------------------------------
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
        #if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g[isave,,] = sigma_gamma_ik
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat
        if(!is.na(match('Ypred', mcmc_params)) || computeDIC) post.Ypred[isave,,] = Yhat + sigma_e*rnorm(length(Y))
        if(!is.na(match('Zhat', mcmc_params)) || computeDIC) post.Zhat[isave,,] = Zhat
        if(!is.na(match('Zpred', mcmc_params)) || computeDIC) post.Zpred[isave,,] = (Yhat + sigma_e*rnorm(length(Y)))^2 * E
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Yhat), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  #if(!is.na(match('sigma_g', mcmc_params))) mcmc_output$sigma_g = post.sigma_g
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_e), m*n),
                            log = TRUE), na.rm=TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}







##### 使用固定的基函数-LR-TPS
#' MCMC Sampling Algorithm for the Count-valued Function-on-Scalars Regression Model
#' with pre-fixed basis functions
#'
#' Runs the MCMC for the function-on-scalars regression model with count responses based on
#' a known basis expansion(splines). Here we assume the factor regression has ICAR errors,
#' as well as some additional default conditions.
#'
#' @param Z the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of FPCs to use; if NULL, select to explain \code{pve}% of variability
#' (ignored for other basis functions)
#' @param pve proportion of variability explained for the FPC basis; only used if \code{K} is NULL
#' @param use_approx_poisson logical; if TRUE, set the size parameter \code{r = 1000} to approximate the Poisson distribution
#' @param sample_r logical; if TRUE, sample the size parameter \code{r},
#' otherwise compute a method-of-moments estimator (unless \code{use_approx_poisson = TRUE})
#' @param sample_sigma logical; if TRUE, sample the variance parameter \code{sigma},
#' otherwise use \code{sigma = 0.1} to include a small amount of small-scale variability
#' @param Offset the \code{T x m} matrix of offsets; if NULL, assume no offsets
#' @param W adjacent matrix
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept)
#' \item "r" (Negative Binomial size parameter)
#' \item "sigma_e" (factor-level error standard deviation)
#' \item "Zhat" (conditional expectation of the count observations, Z)
#' \item "Zpred" (posterior predictive distribution of Z)
#' \item "Zfore" (one-step forecast)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @examples
#' \dontrun{
#' # Fixme
#'}
#'
#' @import KFAS truncdist dfosr
#' @importFrom refund fpca.face
#' @importFrom BayesLogit rpg
# @export
fosr.count_icar_spline = function(Z, tau, X = NULL,  W,
                           use_approx_poisson = FALSE, sample_r = TRUE, sample_sigma = TRUE,
                           Offset = NULL,
                           nsave = 1000, nburn = 1000, nskip = 2,
                           mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                           computeDIC = TRUE){

  # Some options (for now):

  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?

  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100

  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:

 Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Z)))), orthonormalize = TRUE)$Bmat
 #Fmat = Fmat[,1:30]

  # Initialize the FLC coefficients and factors:
  #inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta;  Psi = inits$Psi
     #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #K = ncol(Beta) # to be sure we have the right value


  Beta = tcrossprod(Ztilde, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value


  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)

  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1

  # Impute the missing points:
  #if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
  #                                    size = r,
  #                                    prob = 1 - Pi[na.ind])


  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)


  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))

  # Update the error SD for eta:
  sigma_gamma_tk = matrix(rep(sigma_delta_k, each = T)/sqrt(rep(rowSums(W), K)), nrow = T)

  #---------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term

    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term

    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))

    # Transform for FDLM parameters:
    BtTheta = tcrossprod(Theta - logOffset, t(Fmat))
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------


    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = BtTheta; sigma_tilde = sigma_e

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 )
      gamma_k = gamma_tk[,k]

      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_tk:

    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_tk^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2

    for(t in 1:T){
      ell_gamma_t = ell_gamma_tk1[t,] + (W %*% gamma_tk)[t,]/sigma_delta_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_t * postSD[t,]^2   ##逐个元素相乘 只有一行

      gamma_tk[t, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[t,]))  ##按行更新，迭代新值进入下一行的更新
    }


    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk

    ##update eta
    eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)

    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)

    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(eta_tk*sqrt(rep(rowSums(W),K)), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Sample the corresponding prior variance term(s):

    # Update the error SD for gamma:
    sigma_gamma_tk = rep(sigma_delta_k, each = T)/sqrt( matrix(rep(rowSums(W),K),nrow=T) )
    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,1] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}








##### 基函数-LR-TPS + NIG priors

#' MCMC Sampling Algorithm for the Count-valued Function-on-Scalars Regression Model
#' with pre-fixed basis functions
#'
#' Runs the MCMC for the function-on-scalars regression model with count responses based on
#' a known basis expansion(splines). Here we assume the factor regression has ICAR errors,
#' as well as some additional default conditions.
#'
#' @param Z the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param use_approx_poisson logical; if TRUE, set the size parameter \code{r = 1000} to approximate the Poisson distribution
#' @param sample_r logical; if TRUE, sample the size parameter \code{r},
#' otherwise compute a method-of-moments estimator (unless \code{use_approx_poisson = TRUE})
#' @param sample_sigma logical; if TRUE, sample the variance parameter \code{sigma},
#' otherwise use \code{sigma = 0.1} to include a small amount of small-scale variability
#' @param Offset the \code{T x m} matrix of offsets; if NULL, assume no offsets
#' @param W adjacent matrix
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept)
#' \item "r" (Negative Binomial size parameter)
#' \item "sigma_e" (factor-level error standard deviation)
#' \item "Zhat" (conditional expectation of the count observations, Z)
#' \item "Zpred" (posterior predictive distribution of Z)
#' \item "Zfore" (one-step forecast)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @examples
#' \dontrun{
#' # Fixme
#'}
#'
#' @import KFAS truncdist dfosr
#' @importFrom refund fpca.face
#' @importFrom BayesLogit rpg
# @export
fosr.count_icar_spline_nig = function(Z, tau, X = NULL,  W,
                                  use_approx_poisson = FALSE, sample_r = TRUE, sample_sigma = TRUE,
                                  Offset = NULL,
                                  nsave = 1000, nburn = 1000, nskip = 2,
                                  mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                                  computeDIC = TRUE){

  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100

  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:

  Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Z)))), orthonormalize = TRUE)$Bmat
  #Fmat = Fmat[,1:60]

  # Initialize the FLC coefficients and factors:
  #inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta;  Psi = inits$Psi
  #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #K = ncol(Beta) # to be sure we have the right value


  Beta = tcrossprod(Ztilde, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value


  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)

  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1

  # Impute the missing points:
  #if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
  #                                    size = r,
  #                                    prob = 1 - Pi[na.ind])


  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])


  # SD term for mu_k:
  sigma_mu_k = rep(100, K) # This will stay fixed
  #----------------------------------------------------------------------------

  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)

  # SD term for gamma_tk:
  # Update the error SD for eta:
  sigma_gamma_k = apply(eta_tk, 2, sd)
  sigma_gamma_tk = matrix(rep(apply(eta_tk, 2, sd), each = T)/sqrt(rep(rowSums(W), K)), nrow = T)

  #----------------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
    sigma_omega_pk= matrix(rep(apply(omega, 1, sd), times = K), nrow = p-1)
  }

  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))

    # Transform for FDLM parameters:
    BtTheta = tcrossprod(Theta - logOffset, t(Fmat))
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------


    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = BtTheta; sigma_tilde = sigma_e

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 )
      gamma_k = gamma_tk[,k]

      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_tk:

    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_tk^2,nrow = T))
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2

    for(t in 1:T){
      ell_gamma_t = ell_gamma_tk1[t,] + (W %*% gamma_tk)[t,]/sigma_gamma_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_t * postSD[t,]^2   ##逐个元素相乘 只有一行

      gamma_tk[t, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[t,]))  ##按行更新，迭代新值进入下一行的更新
    }


    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk

    ##update eta
    eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)

    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)

    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # SD term for gamma_tk:

    # Variance part:
    sigma_gamma_k =  apply(eta_tk, 2, function(x){
      1/sqrt(rgamma(n = 1, shape = 0.001 + T/2, rate = 0.001 + sum(x^2 * rowSums(W))/2))
    })


     # Sample the corresponding prior variance term(s):
    # Update the error SD for gamma:
    sigma_gamma_tk = rep(sigma_gamma_k, each = T)/sqrt( matrix(rep(rowSums(W),K),nrow=T) )


    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:

    if(p > 1){
      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      # SD term for omega_pk:
      sigma_omega_pk = matrix(rep(apply(omega, 1, function(x){
        1/sqrt(rgamma(n = 1, shape = 0.001 + K/2, rate = 0.001 + sum(x^2)/2))
      }), times = K), nrow = p-1)
    }


    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,1] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}







######################### NIG priors
#' MCMC Sampling Algorithm for the Count-valued Function-on-Scalars Regression Model
#' with NIG priors
#'
#' Runs the MCMC for the function-on-scalars regression model with count responses based on
#' a reduced-rank expansion. Here we assume the factor regression has ICAR errors,
#' as well as some additional default conditions.
#'
#' @param Z the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param use_approx_poisson logical; if TRUE, set the size parameter \code{r = 1000} to approximate the Poisson distribution
#' @param sample_r logical; if TRUE, sample the size parameter \code{r},
#' otherwise compute a method-of-moments estimator (unless \code{use_approx_poisson = TRUE})
#' @param sample_sigma logical; if TRUE, sample the variance parameter \code{sigma},
#' otherwise use \code{sigma = 0.1} to include a small amount of small-scale variability
#' @param Offset the \code{T x m} matrix of offsets; if NULL, assume no offsets
#' @param W adjacent matrix
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept)
#' \item "r" (Negative Binomial size parameter)
#' \item "sigma_e" (factor-level error standard deviation)
#' \item "Zhat" (conditional expectation of the count observations, Z)
#' \item "Zpred" (posterior predictive distribution of Z)
#' \item "Zfore" (one-step forecast)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @examples
#' \dontrun{
#' # Fixme
#'}
#'
#' @import KFAS truncdist dfosr
#' @importFrom BayesLogit rpg
# @export
fosr.count_icar_nig = function(Z, tau, X = NULL, K = NULL, W,
                           use_approx_poisson = TRUE, sample_r = TRUE, sample_sigma = TRUE,
                           Offset = NULL,
                           nsave = 1000, nburn = 1000, nskip = 2,
                           mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                           computeDIC = TRUE){


  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    #Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100

  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)

  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1

  # Impute the missing points:
  if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                      size = r,
                                      prob = 1 - Pi[na.ind])

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)

  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  sigma_mu_k = rep(100, K) # This will stay fixed
  #----------------------------------------------------------------------------

  #ICAR spatial error eta_ki:  ICAR model(difference)
  eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)

  # SD term for gamma_tk:
  # Update the error SD for eta:
  sigma_gamma_k = apply(eta_tk, 2, sd)
  sigma_gamma_tk = matrix(rep(apply(eta_tk, 2, sd), each = T)/sqrt(rep(rowSums(W), K)), nrow = T)

  #----------------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
    sigma_omega_pk= matrix(rep(apply(omega, 1, sd), times = K), nrow = p-1)
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))

    # Transform for FDLM parameters:
    BtTheta = tcrossprod(t(splineInfo$Bmat), Theta - logOffset)
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtTheta,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_e^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtTheta, Psi); sigma_tilde = sigma_e

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 )   ###修改
      gamma_k = gamma_tk[,k]

      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, (y_tilde_k- gamma_k )/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_tk:

    # postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_tk^2))
    # postMean = matrix((Y_tilde - X%*%alpha_pk)/rep(sigma_tilde^2, times = K))*postSD^2   ##matrix form
    # gamma_tk = matrix(rnorm(n = T*K, mean = postMean, sd = postSD), nrow = T)


    ######
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + 1/sigma_gamma_tk^2)
    ell_gamma_tk1=(Y_tilde - X%*%alpha_pk)/sigma_tilde^2

    for(t in 1:T){
      ell_gamma_t = ell_gamma_tk1[t,] + (W %*% gamma_tk)[t,]/sigma_gamma_k^2  ###ell_gamma_tk的某一行
      postMean = ell_gamma_t * postSD[t,]^2   ##逐个元素相乘 只有一行

      gamma_tk[t, ]=t(rnorm(n=K, mean= postMean,  sd=postSD[t,]))  ##按行更新，迭代新值进入下一行的更新
    }


    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk

    ##update eta
    eta_tk = gamma_tk - (W %*% gamma_tk) * matrix(rep(1/rowSums(W),K), T)

    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)

    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # Variance part:
    # SD term for gamma_tk:

    # Variance part:
    sigma_gamma_k =  apply(eta_tk, 2, function(x){
      1/sqrt(rgamma(n = 1, shape = 0.001 + T/2, rate = 0.001 + sum(x^2 * rowSums(W))/2))
    })


    # Sample the corresponding prior variance term(s):
    # Update the error SD for gamma:
    sigma_gamma_tk = rep(sigma_gamma_k, each = T)/sqrt( matrix(rep(rowSums(W),K),nrow=T) )

    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:

    if(p > 1){
      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      # SD term for omega_pk:
      sigma_omega_pk = matrix(rep(apply(omega, 1, function(x){
        1/sqrt(rgamma(n = 1, shape = 0.001 + K/2, rate = 0.001 + sum(x^2)/2))
      }), times = K), nrow = p-1)
    }

    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}







######################### 不考虑空间效应，没有eta_tk
#' MCMC Sampling Algorithm for the Count-valued Function-on-Scalars Regression Model
#'
#' Runs the MCMC for the function-on-scalars regression model with count responses based on
#' a reduced-rank expansion. Here we assume the factor regression has independent errors,
#' as well as some additional default conditions.
#'
#' @param Z the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param use_approx_poisson logical; if TRUE, set the size parameter \code{r = 1000} to approximate the Poisson distribution
#' @param sample_r logical; if TRUE, sample the size parameter \code{r},
#' otherwise compute a method-of-moments estimator (unless \code{use_approx_poisson = TRUE})
#' @param sample_sigma logical; if TRUE, sample the variance parameter \code{sigma},
#' otherwise use \code{sigma = 0.1} to include a small amount of small-scale variability
#' @param Offset the \code{T x m} matrix of offsets; if NULL, assume no offsets
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept)
#' \item "r" (Negative Binomial size parameter)
#' \item "sigma_e" (factor-level error standard deviation)
#' \item "Zhat" (conditional expectation of the count observations, Z)
#' \item "Zpred" (posterior predictive distribution of Z)
#' \item "Zfore" (one-step forecast)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @examples
#' \dontrun{
#' # Fixme
#'}
#'
#' @import KFAS truncdist dfosr
#' @importFrom BayesLogit rpg
# @export
fosr.count_ind = function(Z, tau, X = NULL, K = NULL,
                           use_approx_poisson = TRUE, sample_r = TRUE, sample_sigma = TRUE,
                           Offset = NULL,
                           nsave = 1000, nburn = 1000, nskip = 2,
                           mcmc_params = list("beta", "fk", "alpha", "Zhat", "Zpred"),
                           computeDIC = TRUE){

  # Some options (for now):
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?
  sample_nu = TRUE # Sample DF parameter, or fix at nu=3?
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  T = nrow(Z); m = ncol(Z); d = ncol(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Check for missing values:
  Zna = Z # The original data, including NAs
  any.missing = any(is.na(Zna));
  if(any.missing) {
    na.ind = which(is.na(Zna), arr.ind = TRUE)
    #Z[na.ind] = mean(Z, na.rm=TRUE) # Initialize at the mean, which should be cleaner
  }
  #----------------------------------------------------------------------------
  # Values for r:
  if(use_approx_poisson){
    r = 1000; # Large value for r approximates the Poisson
    sample_r = FALSE # to be sure
  } else {
    # Method of moments estimator for r:
    r = mean(Z, na.rm=TRUE)^2/(sd(Z, na.rm=TRUE)^2 - mean(Z, na.rm=TRUE))
  }
  # If sampling r, use r~ C+(0, 1/sqrt(e0)) w/ hyperparameter e0:
  if(sample_r)  e0 = 1/100

  # Acting value of Z (for now): simulate from the prior
  Ztilde = log(Z + 1)
  #xi_tj = matrix(NA, nrow = T, ncol = m); not.na.ind = which(!is.na(Z), arr.ind = TRUE)
  #xi_tj[not.na.ind] = rpg(num = nrow(not.na.ind), h = Z[not.na.ind] + r, z = -log(r))
  #Ztilde = (Z - r)/(2*xi_tj) + log(r)
  #----------------------------------------------------------------------------
  # Initialize the main terms:

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Ztilde, tau, K, use_pace = TRUE); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  #inits = fdlm_init_d(Ztilde, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Predicted Mean Values:
  Mu = Theta = tcrossprod(Beta, Fmat)
  Pi = exp(Theta)/(r + exp(Theta))   #Pi = invlogit(Theta)

  # Standard deviation at FDLM level
  if(sample_sigma){
    sigma_e = sd(Ztilde - Mu, na.rm=TRUE); px_sigma_e = 1
  } else sigma_e = 0.1

  # Impute the missing points:
  if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                      size = r,
                                      prob = 1 - Pi[na.ind])

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)

  # Construct the log-offset:
  if(is.null(Offset)){logOffset = 0} else logOffset = log(Offset)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)

    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])

    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_e,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Theta, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # Initialize the corresponding SD term(s):
  xi_gamma_tk = 1/gamma_tk^2; # Precision scale
  nu = 3  # (initial) degrees of freedom


  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K);
  sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))

  # Update the error SD for eta / conditional SD for gamma:
  sigma_gamma_tk = rep(sigma_delta_k, each = T)/sqrt(xi_gamma_tk)

  #---------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term

    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term

    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post.r = array(NA, c(nsave, 1))
  if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi = array(NA, c(nsave, T, m))
  if(!is.na(match('Zhat', mcmc_params))) post.Zhat = array(NA, c(nsave, T, m))
  if(!is.na(match('Zpred', mcmc_params))) post.Zpred = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Step 1: Impute Zt(tau)
    #----------------------------------------------------------------------------
    if(any.missing) Z[na.ind] = rnbinom(n = nrow(na.ind),
                                        size = r,
                                        prob = 1 - Pi[na.ind])
    #----------------------------------------------------------
    # Step 2: Sample r (if desired)
    #----------------------------------------------------------
    if(sample_r){
      expTheta = exp(Theta)
      r = uni.slice(r, g = function(x){
        sum(dnbinom(Z, size = x, prob = x/(x + expTheta), log = TRUE)) - log(1 + x^2*e0)
        #sum(lgamma(Z + x)) - T*m*lgamma(x) + x*sum(log(1-Pi)) - log(1 + x^2*e0)
      }, lower = 0, upper = 1000)
    }
    #----------------------------------------------------------
    # Step 3: PX sample of PG variates and transform to pseudo-Gaussian space
    #----------------------------------------------------------
    xi_tj = rpg(num = T*m, h = Z + r, z = Theta - log(r))
    Ztilde = (Z - r)/(2*xi_tj)  + log(r)
    #----------------------------------------------------------
    # Step 4: Sample Theta
    #----------------------------------------------------------
    postSD = 1/sqrt(xi_tj + 1/sigma_e^2)
    postMean = (matrix(Ztilde)*xi_tj + matrix(Mu)/sigma_e^2)*postSD^2
    Theta = matrix(rnorm(n = T*m, mean = postMean, sd = postSD), nrow = T)
    Pi = exp(Theta)/(r + exp(Theta))

    # Transform for FDLM parameters:
    BtTheta = tcrossprod(t(splineInfo$Bmat), Theta - logOffset)
    #----------------------------------------------------------------------------
    # Step 5: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtTheta,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_e^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 6: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtTheta, Psi); sigma_tilde = sigma_e

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 + sigma_gamma_tk[,k]^2)

      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < T (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, y_tilde_k/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_tk:
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_tk^2))
    postMean = matrix((Y_tilde - X%*%alpha_pk)/rep(sigma_tilde^2, times = K))*postSD^2
    gamma_tk = matrix(rnorm(n = T*K, mean = postMean, sd = postSD), nrow = T)


    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk

    # Update Mu:
    Mu = logOffset + tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 7: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(sample_sigma){
      # Sample the variance:
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + T*m/2,
                              rate = px_sigma_e + sum((Theta - Mu)^2, na.rm=TRUE)/2))
      # Parameter-expanded piece
      px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)

    }
    #----------------------------------------------------------------------------
    # Step 8: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(gamma_tk*sqrt(xi_gamma_tk), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Sample the corresponding prior variance term(s):
    xi_gamma_tk = matrix(rgamma(n = T*K,
                                shape = nu/2 + 1/2,
                                rate = nu/2 + (gamma_tk/rep(sigma_delta_k, each = T))^2/2), nrow = T)
    # Sample degrees of freedom?
    if(sample_nu){
      nu = uni.slice(nu, g = function(nu){
        sum(dgamma(xi_gamma_tk, shape = nu/2, rate = nu/2, log = TRUE)) +
          dunif(nu, min = 2, max = 128, log = TRUE)}, lower = 2, upper = 128)
    }

    # Update the error SD for gamma:
    sigma_gamma_tk = rep(sigma_delta_k, each = T)/sqrt(xi_gamma_tk)
    #----------------------------------------------------------------------------
    # Step 9: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 10: Store the MCMC output:
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('r', mcmc_params)) || computeDIC) post.r[isave,] = r
        if(!is.na(match('Pi', mcmc_params)) || computeDIC) post.Pi[isave,,] = Pi
        if(!is.na(match('Zhat', mcmc_params))) post.Zhat[isave,,] = exp(Mu + sigma_e^2/2)
        if(!is.na(match('Zpred', mcmc_params))) {
          if(use_approx_poisson){
            #post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Theta))
            post.Zpred[isave,,] = rpois(n = T*m, lambda = exp(Mu + sigma_e^2/2))
          } else post.Zpred[isave,,] = rnbinom(n = T*m, size = r, prob = 1 - Pi)
        }
        if(computeDIC) post_loglike[isave] = sum(dnbinom(matrix(Zna), size = r, prob = matrix(1-Pi), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post.r
  if(!is.na(match('Pi', mcmc_params))) mcmc_output$Pi = post.Pi
  if(!is.na(match('Zhat', mcmc_params))) mcmc_output$Zhat = post.Zhat
  if(!is.na(match('Zpred', mcmc_params))) mcmc_output$Zpred = post.Zpred

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnbinom(matrix(Zna),
                              size = colMeans(post.r),
                              prob = matrix(1-colMeans(post.Pi)),
                              log = TRUE), na.rm = TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  return (mcmc_output);
}



