rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)
library(rlist)
library(limSolve)

# load model and dataset
load("inst/logistic/germancredit.RData")

# I replaced the parallel RNG of L'Ecuyer et al (2002) that is used in the previous code
# due to trouble using it on my system. Now it's just set.seed(irep) for iteration irep

# no. of repetitions desired for each processor
nreps <- 30*2^5 # 30 repeats with R = 2^5

# guideline for the choice of k and m values for this application
k <- 330
m <- 10 * k

# stepsize and no. of steps with good contraction and small average compute time
stepsize <- 0.0125
nsteps <- 10

# define HMC kernel and coupled HMC kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define RWMH kernel and coupled RWMH kernel
omega <- 1 / 20 # probability of selecting coupled RWMH
Sigma_std <- 1e-3 # proposal standard deviation of RWMH
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension, gradlogtarget)

# define mixture kernel
mixture_kernel <- function(chain_state, current_pdf, iteration, chain_grad){
  if (runif(1) < omega){
    return(mh$kernel(chain_state, current_pdf, iteration, chain_grad))
  } else {
    return(hmc$kernel(chain_state, current_pdf, iteration, chain_grad))
  }
}

# define coupled mixture kernel
mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration, chain_grad1, chain_grad2){
  if (runif(1) < omega){
    return(mh$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration, chain_grad1, chain_grad2))
  } else {
    return(hmc$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration, chain_grad1, chain_grad2))
  }
}


############################################################################################################
# The following code is somewhat memory-intensive.
# The two loops could easily be combined into one so that the full nreps iterations don't need to be stored

# Part one: Running the coupled HMC chains
############################################################################################################

# pre-allocate
runtimes <- rep(0, nreps)
meetingtime <- rep(0, nreps)
iteration <- rep(0, nreps)
unbiased_estimates <- mcmc_estimates <- matrix(nrow = nreps, ncol = 2 * dimension) # first and second moment

y <- list()
y0 <- list()
weights <- list()
mcmc_ind <- list()
chain <- list()
X_for_ZV <- list()
X_for_ZV0 <- list()


for(irep in 1:nreps){
  set.seed(irep)

  tic()
    estimation_output <- unbiased_estimator(logtarget, mixture_kernel, mixture_coupled_kernel,
                                            rinit, h = function(x) c(x, x^2), k = k, m = m, gradlogtarget = gradlogtarget)
  timing <- toc()
  runtime <- timing$toc - timing$tic
  runtimes[irep] <- runtime
  meetingtime[irep] <- estimation_output$meetingtime
  iteration[irep] <- estimation_output$iteration
  unbiased_estimates[irep, ] <- estimation_output$uestimator
  mcmc_estimates[irep, ] <- estimation_output$mcmcestimator
  y[[irep]] <- estimation_output$y
  weights[[irep]] <- estimation_output$weights
  chain[[irep]] <- estimation_output$chain
  mcmc_ind[[irep]] <- estimation_output$mcmc_ind

  # First order polynomial ZV-CV design matrix. The factor of 0.5 is inconsequential but is used to match the Mira et al. (2013) ZV-CV paper.
  # Higher order polynomials could be obtained using the getX function from the ZVCV CRAN package.
  X_for_ZV[[irep]] <- cbind(rep(1,length(weights[[irep]])),-0.5*estimation_output$gradients)

  cat("Repetition:", irep, "/", nreps, "\n")
}

filename <- paste("output.hmc.germancredit.RData", sep = "")
save(nreps, k, m, stepsize, nsteps, runtimes, meetingtime, iteration, unbiased_estimates, mcmc_estimates,
     y, y0, weights, chain, X_for_ZV, X_for_ZV0, mcmc_ind, repeat_index,
     file = filename, safe = F)



############################################################################################################
# Part two: Performing the variance reduction
############################################################################################################

R <- 2^5 # The code is set up for this being divisible by 2 (adjustments to non-even R is simple)
total_num <- floor(nreps/R)

tt <- R/2 - 1
inds <- seq(1,nreps-R+1,by=R)

unbiased_combined <- matrix(NaN,nrow=total_num,ncol=2*dimension)
unbiased_ub_combined <- matrix(NaN,nrow=total_num,ncol=2*dimension)
unbiased_ev_combined <- matrix(NaN,nrow=total_num,ncol=2*dimension)
for (k in 1:length(inds)){
  unbiased_combined[k,] <- colMeans(unbiased_estimates[(inds[k]):(inds[k]+2*tt+1),])
}

# Empirical variance of the estimators
optim_fun <- function(beta,est_plain,X,y,weights){
  est_new <- est_plain
  for (i in 1:length(est_plain)){
    est_new[i] <- est_new[i] - weights[[i]]%*%X[[i]][,2:(length(beta)+1)]%*%matrix(beta,nrow=length(beta))
  }
  return(mean((est_new - mean(est_new))^2))
}

# Gradient of the empirical variance of the estimators wrt the parameters beta
optim_grad_fun <- function(beta,est_plain,X,y,weights){
  R <- length(est_plain)
  est_new <- est_plain
  grad_temp <- matrix(NaN,nrow=R,ncol=length(beta))
  for (i in 1:R){
    est_new[i] <- est_new[i] - weights[[i]]%*%X[[i]][,2:(length(beta)+1)]%*%matrix(beta,nrow=length(beta))
    grad_temp[i,] <- - weights[[i]]%*%X[[i]][,2:(length(beta)+1)]
  }
  return( 2*colMeans((grad_temp - matrix(rep(colMeans(grad_temp),R),nrow=R,byrow=TRUE))*(matrix(rep(est_new,length(beta)),nrow=R) - matrix(rep(mean(est_new),length(beta)*R),nrow=R))) )
}

for (k in 1:length(inds)){

  i <- inds[k]

  # Combining samples together. a is for the first group and b is for the second
  # Group a estimates the function with the first R/2 samples and evaluates the estimator with the second half, and vice versa for b.

  X_a <- list.rbind(X_for_ZV[i:(i+tt)])
  X_b <- list.rbind(X_for_ZV[(i+tt+1):(i+2*tt+1)])
  y_a <- list.rbind(y[i:(i+tt)])
  y_b <- list.rbind(y[(i+tt+1):(i+2*tt+1)])
  weights_a <- unlist(weights[i:(i+tt)])
  weights_b <- unlist(weights[(i+tt+1):(i+2*tt+1)])

  # Least squares estimator
  # This is estimator (ii) of our comment to the Jacob et al JRSS B paper and it's based on minimising the upper bound.

  fit_a <- lm(y_a ~ X_a - 1)
  fit_b <- lm(y_b ~ X_b - 1)

  beta_a_ub <- coef(fit_a)
  beta_b_ub <- coef(fit_b)

  estim_a_ub <- 1/(tt+1)*colSums(weights_a%*%X_a[,2:(dimension+1)]%*%beta_b_ub[2:(dimension+1),])
  estim_b_ub <- 1/(tt+1)*colSums(weights_b%*%X_b[,2:(dimension+1)]%*%beta_a_ub[2:(dimension+1),])

  unbiased_ub_combined[k,] <- 1/(2*tt+2)*colSums(unbiased_estimates[i:(i+2*tt+1),]) - 0.5*(estim_a_ub + estim_b_ub)

  # Minimising the empirical variance (estimator (i) of our comment to the Jacob et al JRSS B paper), staring at the solution based on the upper bound
  # Unfortunately this optimisation seems to be highly multimodal and its performance is strongly dependent on the starting point
  beta_a_ev <- beta_b_ev <- matrix(NaN,nrow=dimension,2*dimension)
  for (j in 1:(2*dimension)){
    init <- beta_a_ub[2:(dimension+1),j]
    beta_a_ev[,j] <- optim(init, optim_fun, gr = optim_grad_fun, est_plain = unbiased_estimates[i:(i+tt),j],X=X_for_ZV[i:(i+tt)],y=y[i:(i+tt)],weights=weights[i:(i+tt)],method='BFGS',control=list(maxit=100))$par

    init <- beta_b_ub[2:(dimension+1),j]
    beta_b_ev[,j] <- optim(init, optim_fun, gr = optim_grad_fun, est_plain = unbiased_estimates[(i+tt+1):(i+2*tt+1),j],X=X_for_ZV[(i+tt+1):(i+2*tt+1)],y=y[(i+tt+1):(i+2*tt+1)],weights=weights[(i+tt+1):(i+2*tt+1)],method='BFGS',control=list(maxit=100))$par

    print(sprintf("--Dimension j index is %d/%d",j,2*dimension))
  }

  estim_a_ev <- 1/(tt+1)*colSums(weights_a%*%X_a[,2:(dimension+1)]%*%beta_b_ev)
  estim_b_ev <- 1/(tt+1)*colSums(weights_b%*%X_b[,2:(dimension+1)]%*%beta_a_ev)

  unbiased_ev_combined[k,] <- 1/(2*tt+2)*colSums(unbiased_estimates[i:(i+2*tt+1),]) - 0.5*(estim_a_ev + estim_b_ev)

  print(sprintf("- Replicate k index is %d",k))
}

filename <- paste("output.hmc.germancredit.cv.RData", sep = "")
save(nreps, k, m, stepsize, nsteps, unbiased_estimates, mcmc_estimates, unbiased_combined, unbiased_ub_combined, unbiased_ev_combined, 
     file = filename, safe = F)
