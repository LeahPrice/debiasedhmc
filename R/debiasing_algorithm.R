# This file consists of the general-purpose functions coupled_chains, H_bar and unbiased_estimator,
# which implement our debiased MCMC algorithm for general kernels and test functions h(.)

# Run coupled chains until max(tau, m) where tau is the meeting time and m specified by user
#'@rdname coupled_chains
#'@title Coupled MCMC chains
#'@description Sample two MCMC chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, m), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'  See \code{\link{get_hmc_kernel}}
#' for an example of function returning the appropriate kernels.
#'@param logtarget function evaluating the log target density
#'@param single_kernel function taking a state (in a vector), its log density and an iteration, and returning
#' a list with a key named \code{chain_state} containing the next state and its log density \code{current_pdf}.
#'@param coupled_kernel function taking two states (in two vectors), their log densities and an iteration,
#' and returning a list with keys \code{chain_state1}, \code{chain_state2}, \code{current_pdf1} and \code{current_pdf2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param m number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{m},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@param gradlogtarget function to compute gradient of target log-density
#'@export
coupled_chains <- function(logtarget, single_kernel, coupled_kernel, rinit, m = 1, max_iterations = Inf, preallocate = 10, gradlogtarget = NULL){
  # keep track of stuff
  # initialize
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  current_pdf1 <- logtarget(chain_state1)
  current_pdf2 <- logtarget(chain_state2)
  if (!is.null(gradlogtarget)){
    chain_grad1 <- gradlogtarget(chain_state1)
    chain_grad2 <- gradlogtarget(chain_state2)
  } else{
    chain_grad1 <- chain_grad2 <- NULL
  }

  # pre-allocate
  p <- length(chain_state1)
  samples1 <- matrix(nrow = m+preallocate+1, ncol = p)
  samples2 <- matrix(nrow = m+preallocate, ncol = p)
  nrowsamples1 <- m+preallocate+1
  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  current_nsamples1 <- 1
  iter <- 1

  # move first chain
  res_single_kernel <- single_kernel(chain_state1, current_pdf1, iter)
  chain_state1 <- res_single_kernel$chain_state
  current_pdf1 <- res_single_kernel$current_pdf
  chain_grad1 <- res_single_kernel$chain_grad
  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1

  # iterate
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      # only need to use single kernel after meeting
      res_single_kernel <- single_kernel(chain_state1, current_pdf1, iter, chain_grad1)
      chain_state1 <- res_single_kernel$chain_state
      current_pdf1 <- res_single_kernel$current_pdf
      chain_grad1 <- res_single_kernel$chain_grad
      chain_state2 <- chain_state1
      current_pdf2 <- current_pdf1
      chain_grad2 <- chain_grad1
    } else {
      # use coupled kernel
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iter, chain_grad1, chain_grad2)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf1 <- res_coupled_kernel$current_pdf1
      current_pdf2 <- res_coupled_kernel$current_pdf2
      chain_grad1 <- res_coupled_kernel$chain_grad1
      chain_grad2 <- res_coupled_kernel$chain_grad2

      # check if meeting happens
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }

    # store coupled chains
    if ((current_nsamples1+1) > nrowsamples1){
      new_rows <- nrowsamples1 - 1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    current_nsamples1 <- current_nsamples1 + 1

    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }

  # drop redundant entries
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

# Compute unbiased estimator using output of coupled_chains
#'@rdname H_bar
#'@title Compute unbiased estimators from coupled chains
#'@description Compute the proposed unbiased estimators given the output
#'of the 'coupled_chains' function. The integral of interest is that of the function h,
#'which can be multivariate. The choice of k and m must be such that m is at most the choice
#'made when running coupled_chains, and k must be less than m.
#'@export
H_bar <- function(c_chains, h = function(x) x, k = 0, m = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (m > maxiter){
    print("error: m has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(m+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  if (c_chains$meetingtime <= k + 1){
    # nothing else to add
  } else {
    deltas <- matrix(0, nrow = maxiter - k + 1, ncol = p)
    deltas_term <- rep(0, p)
    for (t in k:min(maxiter-1, c_chains$meetingtime-1)){ # t is as in the report, where the chains start at t=0
      coefficient <- min(t - k + 1, m - k + 1)
      delta_tp1 <- h(c_chains$samples1[t + 1 + 1,]) - h(c_chains$samples2[t+1,]) # the +1's are because R starts indexing at 1
      deltas_term <- deltas_term + coefficient * delta_tp1
    }
    H_bar <- H_bar + deltas_term
  }
  return(H_bar / (m - k + 1))
}

# Run coupled chains until max(tau, m), where tau is the meeting time and m specified by user,
# and compute H_bar estimator simultaneously, without storing the chains
#'@rdname unbiased_estimator
#'@title Unbiased estimator
#'@description Sample two MCMC chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, m), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'  See \code{\link{get_hmc_kernel}}
#' for an example of function returning the appropriate kernels.
#'
#'@param logtarget function evaluating the log target density
#'@param single_kernel function taking a state (in a vector), its log density and an iteration, and returning
#' a list with a key named \code{chain_state} containing the next state and its log density \code{current_pdf}.
#'@param coupled_kernel function taking two states (in two vectors), their log densities and an iteration,
#' and returning a list with keys \code{chain_state1}, \code{chain_state2}, \code{current_pdf1} and \code{current_pdf2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param h test function (possibly vector-valued)
#'@param k burn-in parameter (default to 0)
#'@param m time average parameter (will be proportional to the computing cost if meeting occurs before \code{m},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param gradlogtarget function to compute gradient of target log-density
#'@export
unbiased_estimator <- function(logtarget, single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf, gradlogtarget = NULL){
  # initialize
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  current_pdf1 <- logtarget(chain_state1)
  current_pdf2 <- logtarget(chain_state2)
  if (!is.null(gradlogtarget)){
    chain_grad1 <- gradlogtarget(chain_state1)
    chain_grad2 <- gradlogtarget(chain_state2)
  } else{
    chain_grad1 <- chain_grad2 <- NULL
  }

  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  samples0 <- chain_state1
  gradients0 <- chain_grad1
  fun_evals0 <- mcmcestimator

  samples <- NULL
  gradients <- NULL
  fun_evals <- NULL
  weights <- NULL
  mcmc_ind <- NULL
  chain <- NULL
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  } else{
    samples <- rbind(samples,chain_state1)
    gradients <- rbind(gradients,chain_grad1)
    fun_evals <- rbind(fun_evals,h(chain_state1))
    weights <- c(weights,1 / (m - k + 1))
    mcmc_ind <- c(mcmc_ind,TRUE)
    chain <- c(chain,"x")
  }

  # move first chain
  iter <- 1
  res_single_kernel <- single_kernel(chain_state1, current_pdf1, iter, chain_grad1)
  chain_state1 <- res_single_kernel$chain_state
  current_pdf1 <- res_single_kernel$current_pdf
  chain_grad1 <- res_single_kernel$chain_grad

  # correction term computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  if (k == 0){
    correction <- correction + min(1, (0 - k + 1)/(m - k + 1) )  * ( h(chain_state1) - h(chain_state2) )
    samples <- rbind(samples,chain_state2)
    gradients <- rbind(gradients,chain_grad2)
    fun_evals <- rbind(fun_evals,h(chain_state2))
    weights <- c(weights, - min(1, (0 - k + 1)/(m - k + 1) ))
    mcmc_ind <- c(mcmc_ind,FALSE)
    chain <- c(chain,"y")
    samples <- rbind(samples,chain_state1)
    gradients <- rbind(gradients,chain_grad1)
    fun_evals <- rbind(fun_evals,h(chain_state1))
    weights <- c(weights,min(1, (0 - k + 1)/(m - k + 1) ))
    mcmc_ind <- c(mcmc_ind,FALSE)
    chain <- c(chain,"x")
  }



  # accumulate mcmc estimator
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
    if (k==0){
      weights[length(weights)] <- weights[length(weights)] + 1 / (m - k + 1)
      mcmc_ind[length(mcmc_ind)] <- TRUE
    }
  }

  # iterate
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      # only need to use single kernel after meeting
      res_single_kernel <- single_kernel(chain_state1, current_pdf1, iter, chain_grad1)
      chain_state1 <- res_single_kernel$chain_state
      current_pdf1 <- res_single_kernel$current_pdf
      chain_grad1 <- res_single_kernel$chain_grad
      chain_state2 <- chain_state1
      current_pdf2 <- current_pdf1
      chain_grad2 <- chain_grad1

      # accumulate mcmc estimator
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
        samples <- rbind(samples,chain_state1)
        gradients <- rbind(gradients,chain_grad1)
        fun_evals <- rbind(fun_evals,h(chain_state1))
        weights <- c(weights,1 / (m - k + 1))
        mcmc_ind <- c(mcmc_ind,TRUE)
        chain <- c(chain,"x")
      }

    } else {
      # use coupled kernel
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iter, chain_grad1, chain_grad2)
      chain_state1 <- res_coupled_kernel$chain_state1
      current_pdf1 <- res_coupled_kernel$current_pdf1
      chain_grad1 <- res_coupled_kernel$chain_grad1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf2 <- res_coupled_kernel$current_pdf2
      chain_grad2 <- res_coupled_kernel$chain_grad2

      # check if meeting happens
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }

      if (k <= iter){
        # accumulate mcmc estimator
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
          samples <- rbind(samples,chain_state1)
          gradients <- rbind(gradients,chain_grad1)
          fun_evals <- rbind(fun_evals,h(chain_state1))
          weights <- c(weights,1 / (m - k + 1))
          mcmc_ind <- c(mcmc_ind,TRUE)
          chain <- c(chain,"x")
        }

        # accumulate correction term
        correction <- correction + min(1, (iter-1 - k + 1)/(m - k + 1) ) * ( h(chain_state1) - h(chain_state2) )
        if (iter<=m){
          weights[length(weights)] <- weights[length(weights)] + min(1, (iter-1 - k + 1)/(m - k + 1) )
        } else{
          samples <- rbind(samples,chain_state1)
          gradients <- rbind(gradients,chain_grad1)
          fun_evals <- rbind(fun_evals,h(chain_state1))
          weights <- c(weights,min(1, (iter-1 - k + 1)/(m - k + 1) ))
          mcmc_ind <- c(mcmc_ind,FALSE)
          chain <- c(chain,"x")
        }
        samples <- rbind(samples,chain_state2)
        gradients <- rbind(gradients,chain_grad2)
        fun_evals <- rbind(fun_evals,h(chain_state2))
        weights <- c(weights, - min(1, (iter-1 - k + 1)/(m - k + 1) ))
        mcmc_ind <- c(mcmc_ind,FALSE)
        chain <- c(chain,"y")
      }

    }

    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }

  # compute unbiased estimator
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction

  # mcmcestimator_reproduced <- colMeans(matrix(fun_evals[mcmc_ind,],nrow=sum(mcmc_ind)))
  # uestimator_reproduced <- t(weights%*%fun_evals)

  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              samples = samples, gradients = gradients, y = fun_evals, weights = weights, chain = chain, mcmc_ind = mcmc_ind,
              samples0 = samples0, gradients0 = gradients0, y0 = fun_evals0,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}
