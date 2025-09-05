


# required libraries and functions
require(VGAM)

# logsumexp avoids underflow 
logsumexp = function(x){
  i = which.max(x); 
  res = x[i] + log1p(sum(exp(x[-i] - x[i]))); 
  if (is.nan(res)) res=-Inf; return(res)
}

# this file contains the following functions : em_ft_clust, em_initialization and ft_clust functions
# em_ft_clust and em_initialization functions share same input values such that they are presented together. 

#' em_ft_clust function is the implementation of the EM algorithm for FT-clustering
#' em_initialization function returns data driven initial values for em_ft_clust. It is the implementation of an EM algorithm over a model solely composed of non-neutral groups (s.t. the frequency trajectory of the neutral group has sae distribution as non-neutral groups).
#' @param data a mandatory list containing 
#'  data$x: n.nT matrix of mutation counts where n is the number of mutations and nT is the number of time points. The character vector of mutation names as row names is optional. The character vector of dates as column names is mandatory.
#'  data$d: n.nT matrix of read depths at related positions 
#' @param K single value for the number of non-neutral groups (default: NULL)
#' @param eta n.(K+1) matrix of posterior group assignment where the first column stands for the neutral group (default: NULL)
#' Either eta or K must be provided. If eta is provided, it is used as initial values. If eta is NULL, K is used for a random initialization of eta. 
#' @param niter single value for the maximal number of iterations (default: 2000)
#' @param tol single value for the tolerance on parameters for the algorithm to break, converge (default: 1e-3)
#' @param verbose boolean vector. If TRUE, the log-likelihood and parameter estimates are provided at each iteration (default : TRUE for em_ft_clust and FALSE for em_initialization)
em_ft_clust <- function(data, K = NULL, eta = NULL, niter = 2000, tol = 1e-3, verbose = FALSE) {
  
  if(is.null(eta)) {
    # random initialization
    eta = matrix(NA,n,K+1)
    for (i in 1:n) {tmp = runif(K+1); eta[i,] = tmp / (sum(tmp))}
  } else { 
    # otherwise retrieve K if initial eta is provided
    K = ncol(eta)-1 
  }
  
  x = data$x
  d = data$d
  
  # compute vector of time differences 
  time = as.numeric(julian(as.Date(colnames(x)), origin = as.Date(colnames(x)[1])))
  
  n = nrow(x) # number of mutations 
  nT = length(time) # number of time points
  
  if (K==0) { # no group under selection
    # alpha and beta estimates 
    fit = vglm(cbind(apply(x,1,sum), apply(d-x, 1, sum)) ~ 1,family = betabinomialff())
    alpha = Coef(fit)[1] # careful Coef and not coef function
    beta = Coef(fit)[2] # careful Coef and not coef function
    
    # compute loglik
    loglik = 0
    for (i in 1:n) {
      res = sum(lchoose(d[i,],x[i,]))
      for (a in 0:(sum(x[i,])-1)) res = res + log(alpha+a)
      for (a in 0:(sum(d[i,]-x[i,])-1)) res = res + log(beta+a)
      for (a in 0:(sum(d[i,])-1)) res = res - log(alpha + beta + a)
      loglik = loglik + res
    }
    
    return(list(loglik = loglik, alpha = alpha, beta = beta))
    
  } else { # at least one group under positive or negative selection
    
    # theta.tol introduced to check convergence
    theta.tol = rep(0, 2*K+2)
    
    for (iter in 1:niter) {
      
      # M step
      ## update pi (group proportion)
      pi = apply(eta,2,sum) / n
      
      ## update alpha and beta (intercept for the neutral group)
      weights = eta[,1];
      weights[weights==0]=1e-300 # vglm does not support null weights
      fit = vglm(cbind(apply(x,1,sum), apply(d-x, 1, sum)) ~ 1,family = betabinomialff(), weights = weights)
      alpha = Coef(fit)[1] # careful Coef and not coef function
      beta = Coef(fit)[2] # careful Coef and not coef function

      ## update mu and s (intercept and selection coefficient for non-neutral groups)
      df = data.frame(x = rep(c(t(x)), K), d = rep(c(t(d)), K), time = rep(time, n*K), z = as.factor(rep(1:K, each = n*nT)))
      weights = NULL
      for (k in 1:K) weights = c(weights, rep(eta[,k+1], each = nT))
      if(K>1) {
        fit = glm(cbind(x, d - x) ~ 0 + z + time:z, "binomial", df, weights)
      } else {
        fit = glm(cbind(x, d - x) ~ 1 + time, "binomial", df, weights)
      }
      mu = coefficients(fit)[1:K]
      s = coefficients(fit)[-(1:K)]

      # E step
      eta = matrix(NA, n, K+1)
      for (i in 1:n) {
        ## neutral group
        res = sum(lchoose(d[i,], x[i,]))
        for (a in 0:(sum(x[i,])-1)) res = res + log(alpha + a)
        for (a in 0:(sum(d[i,]-x[i,])-1)) res = res + log(beta + a)
        for (a in 0:(sum(d[i,])-1)) res = res - log(alpha + beta + a)
        eta[i,1] = res + log(pi[1])

        ## non-neutral groups
        for (k in 1:K) {
          prob = exp(mu[k]+s[k]*time)
          prob = prob/(1+prob)
          eta[i,k+1] = sum(dbinom(x[i,], d[i,], prob, log=TRUE)) + log(pi[k+1])
        }
      }
      aux = apply(eta, 1, logsumexp)
      loglik = sum(aux) # log-likelihood
      eta = exp(eta - aux)

      if(verbose) {
        cat("iter = ", iter, "loglik = ", loglik, "pi =", pi, "mu = ", mu, "s = ", s, "alpha =", alpha, "beta = ", beta, "\n")
      }
  
      #check convergence
      if(max(abs(theta.tol - c(mu, s, alpha, beta)) / abs(c(mu, s, alpha, beta))) < tol) break;
      theta.tol = c(mu,s, alpha, beta) 
    }
    
    rownames(eta) = rownames(x)
    size.smallest.group = min(table(apply(eta,1,which.max)))
      
    return(list(loglik = loglik, pi = pi, mu = mu, s = s, alpha = alpha, beta = beta, eta = eta, size.smallest.group = size.smallest.group)) 
  }
}
em_initialization <- function(data, K, eta = NULL, niter = 2000, tol = 1e-3, verbose = FALSE) {
  
  if (!all(dim(data$x) == dim(data$d))) stop("inconsistent data$x and data$d dimensions")
  if (all(is.null(K) & is.null(eta))) stop("either eta and/or K must be provided")
  if (all(!is.null(K) & !is.null(eta))) {
    if (ncol(eta) != K+1) stop("inconsistent dimension of eta and/or K")
  }
  
  x = data$x
  d = data$d
  
  # compute vector of time differences 
  time = as.numeric(julian(as.Date(colnames(x)), origin = as.Date(colnames(x)[1])))
  
  n = nrow(x) # number of mutations
  nT = ncol(x) # number of time points
  
  # initialization
  if(is.null(eta)) {
    # if eta is NULL, random initialization of posterior group assignment
    eta = matrix(NA, n, K+1)
    for (i in 1:n) {
      tmp = runif(K+1); 
      eta[i,] = tmp / (sum(tmp))
    }
  } else { 
    # if initial eta is provided, retrieve K 
    K = ncol(eta) - 1 
  }
  
  # theta.tol introduced to check algorithm convergence
  theta.tol = rep(0, 3*(K+1))
  
  for (iter in 1:niter) {
    
    # M step
    pi = apply(eta,2,sum) / n
    df = data.frame(x = rep(c(t(x)), K+1), d = rep(c(t(d)), K+1), time = rep(time, n*(K+1)), z = as.factor(rep(0:K, each = n*nT)))
    weights = NULL
    for (k in 0:K) weights = c(weights, rep(eta[,k+1], each = nT))
    if (K > 0) {
      fit = glm(cbind(x, d - x) ~ 0 + z + time:z, "binomial", df, weights)
    } else {
      fit = glm(cbind(x, d - x) ~ 1 + time, "binomial", df, weights)
    }
    tmp = coefficients(fit)
    mu = tmp[1:(K+1)]; 
    s = tmp[(K+1)+(1:(K+1))]
    
    # E step
    eta = matrix(0, nrow = n, ncol = K+1)
    for (i in 1:n) for (k in 0:K) {
      prob = exp(mu[k+1] + time * s[k+1]) 
      prob = prob / (1 + prob)
      prob[prob==1.0] = 1-1e-10; prob[prob==0.0] = 1e-10
      eta[i,k+1] = sum(dbinom(x[i,], d[i,], prob, log = TRUE)) + log(pi[k+1])
    }
    aux = apply(eta, 1, logsumexp)
    loglik = sum(aux)
    eta = exp(eta - aux)
    
    if(verbose) {
      cat("loglik=",loglik, "theta.tol=", theta.tol,"\n")
    }
    
    # check convergence
    if(max(abs(theta.tol - c(pi,mu, s)) / abs(c(pi,mu, s))) < tol) break;
    theta.tol = c(pi, mu, s)
    
  }
  return(list(loglik = loglik, pi = pi, mu = mu, s = s, eta = eta))
}


#' ft_clust function combines em_ft_clust and em_initialization for a chosen K
#' @param n.pre.init vector of size 2 for the number of random initial eta for em_initialization performed over a model composed of a total of K+2 (first slot) and K+3 (second slot) groups (default: c(5,5)).
#' @param n.pre.it vector of size 2 for the related number of iterations to perform for em_initialization (default: c(5,5)). 
#' @param n.init vector of size 2 for the number of initial eta for em_ft_clust computed with em_initialization performed over a model composed of a total of K+2 (first slot) and K+3 (second slot) groups (default: c(10,10)). 
#' @param n.it vector of size 2 for the related number of iterations to perform for em_ft_clust during initialization step (default: c(5,5)).
#' @param check.initialization boolean vector. If TRUE, the vector of log-likelihood reached after the n.it iterations of em_ft_clust over the n.init values of eta is reported (default: TRUE)
#' @param verbose boolean vector. If TRUE, the log-likelihood (respectively the log-likelihood and parameter estimates) are porvided at each iteration of em.initilization (respectively em_ft_clust)

#' @return a list of the log-likelihood, parameter estimates, posterior group assignment and check.init vector if check.initialization is TRUE

ft_clust <- function(data, K,
                     n.pre.init = c(5,5), n.pre.it = c(5,5),
                     n.init = c(10,10), n.it = c(5,5), 
                     check.initialization = TRUE, 
                     niter = 2000,
                     tol = 1e-3,
                     verbose = TRUE) {
  
  check.init = NULL
  
  loglik = -Inf # initial log-likelihood 
  
  # propose initial values computed with em_initialization over a model composed of a total of K+2 and K+3 groups
  for (k in which(n.init!=0)) {
    if(verbose) cat("k =",k,"\n")
    i = 1
    while (i <= n.init[k]) {
      
      # first step : compute initial eta for em_ft_clust 
        # perform n.pre.init times n.pre.it iterations over em_initialization and keep eta associated to the best log-likelihood
      pre.loglik = -Inf 
      j=1
      while (j <= n.pre.init[k]) {
        tmp = try(em_initialization(data, K + k, eta = NULL, niter = n.pre.it[k], verbose = FALSE))
        if (class(tmp)!= "try-error") {
          if(verbose) {cat("rep.pre.init =",j, "loglik =", tmp$loglik,"\n")}
          if(tmp$loglik > pre.loglik) {pre.loglik = tmp$loglik; pre.init = tmp}
          j = j+1
        }
      }
      
      # 2nd step : Perform n.it iterations of em_ft_clust over the initial eta selected
      # extract initial neutral group by fusing the k+1 groups of lowest absolute value for s
      idx = sort.int(abs(pre.init$s), index.return = TRUE)$ix
      tmp = try(em_ft_clust(data, 
                            K, 
                            eta = cbind(apply(pre.init$eta[,-tail(idx,K)], 1, sum), pre.init$eta[,tail(idx,K)]), 
                            niter = n.it[k],
                            tol = tol,
                            verbose = verbose))
      # store associated loglikelihood of check.initialization is TRUE
      if (check.initialization) {
        check.init = c(check.init, ifelse(class(tmp)!="try-error", tmp$loglik, "error"))
      }
      
      # keep 'init' as eta if it is  associated to the best loglik
      if(class(tmp)!="try-error") {
        if(tmp$loglik > loglik) {
          loglik = tmp$loglik; init = tmp
        }
      }
      # test subsequent value of initial eta
      i = i + 1
    }
    # test next model (composed of increasing total number of groups)
    k = k + 1
  }
  
  # run algorithm until convergence with init associated to best log-likelihood
  res = try(em_ft_clust(data, K, eta = init$eta, niter = niter, tol = 1e-5, verbose = verbose))
  if(check.initialization) {
    names(check.init) = c(rep("K+2", n.init[1]), rep("K+3", n.init[2]))
    res$check.initialization = check.init
  }
  
  return(res = res)
}



