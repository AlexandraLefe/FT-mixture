
# this file contains the following functions : generate_simu

#' generate_simu generates a simulated dataset with the same model as ft_clust
#' @param n single value for the number of mutations
#' @param time vector of the size of the number of time points for the time differences (usually in days) between sampling date and time origin. First slot is zero
#' @param date.origin the date of origin, given is character in the format "2019-01-01". If NULL (default value) the date "2019-01-01" is used. 
#' @param param a list of parameters composed of 
#'  pi vector of size K+1 for group proportions, where K is the number of non-neutral groups
#'  mu vector of size K for the logit of the frequency trajectory of non-neutral groups
#'  s vector of size K for selection coefficients of non-neutral groups
#'  alpha and beta : single values for the parameters of the betabinomial distribution of the intercept of the neutral group.
#' @param lambda single value for the mean of the poisson distribution for read depths
#' 
#' @return a list of a dataset (list of a matrix mutation count and a matrix of read depth), a vector z of group assignment and a list of parameters
generate_simu <- function(n, time, date.origin = NULL, param, lambda) {
  
  pi = param$pi
  mu = param$mu
  s = param$s
  alpha = param$alpha
  beta = param$beta
  
  
  # number of time points 
  nT = length(time)
  # number of non-neutral groups 
  K = length(pi) - 1
  
  # group assignment
  z = NULL; while(length(unique(z)) != K+1) {z = sample(0:K, n, replace = TRUE, prob = pi)}
  
  # matrix of read depth
  d = matrix(rpois(n * nT, lambda = lambda), n, nT)
  
  # matrix of mutation counts
  x = matrix(NA,n,nT)
  for (i in 1:n) {
    k = z[i]
    if (k==0) { 
      x[i,] = rbinom(nT,d[i,],rbeta(1,alpha,beta)) 
    } else {
      x[i,] = rbinom(nT,d[i,],exp(mu[k] + s[k]*time) / (1+exp(mu[k] + s[k]*time)))
    }
  }
  
  # add names
  colnames(x) = colnames(d) = as.character(as.Date (ifelse(is.null(date.origin), "2019-01-01", date.origin)) + time)
  rownames(x) = rownames(d) = paste0("mut.",1:n)
  data = list(x = x, d = d)

  return(list(data = list(x=x, d=d), z = z))
}

# # uncomment to plot raw trajectories for verification
# plot(0, type = "n", xlim = range(time), ylim = c(0,1), xlab = "days", ylab = "frequency")
# for (i in 1:n) {
#   lines(time, (simu.dataset$data$x/simu$data$d)[i,], col = simu$z[i]+1)
# }


require(gtools)
#' handle_label_switching uses the mean square error between true (used) for simulations parameters mu and s and parameter estimates to handle label switching in simulation studies. 
#' @param res an output of ft.clust function over a simulated dataset
#' @param param the parameter used for simulation
handle_label_switching <- function(res,param) {
  K = length(param$pi) - 1
  perm = permutations(K, K, repeats.allowed = FALSE) 
  idx = perm[which.min(apply((matrix(res$mu[perm], K, K, byrow = TRUE) - 
                                matrix(param$mu, K, K, byrow = TRUE))^2, 1, sum) +
                         apply((matrix(res$s[perm]*max(time), K, K, byrow = TRUE) -
                                  matrix(param$s*max(time), K, K, byrow = TRUE))^2, 1, sum)),] 
  res$pi = res$pi[c(1,idx+1)];
  res$mu = res$mu[idx];
  res$s = res$s[idx]; 
  res$eta = res$eta[,c(1,idx+1)]
  
  return(list(res=res, param=param))
}



#' simu_hidden_RW generates a simulated dataset with the hidden random walk model 
#' @param n single value for the number of mutations
#' @param time vector of the size of the number of time points for the time differences (usually in days) between sampling date and time origin. First slot is zero
#' @param date.origin the date of origin, given is character in the format "2019-01-01". If NULL (default value) the date "2019-01-01" is used. 
#' @param param a list of parameters composed of 
#'  pi: vector of size K+1 for group proportions, where K is the number of non-neutral groups
#'  mu: vector of size K for the logit of the frequency trajectory of non-neutral groups
#'  s: vector of size K for selection coefficients of non-neutral groups
#'  alpha and beta: single values for the parameters of the betabinomial distribution of the intercept of the neutral group.
#'  sd0: single value for the standard deviation for the intercept of non-neutral groups
#'  sd: single value for the standard deviation for gaussian random walk (common for all groups). 
#' @param lambda single value for the mean of the poisson distribution for read depths
#' 
#' @return a list of a dataset (list of a matrix mutation count and a matrix of read depth), a vector z of group assignment and a list of parameters
#' 
#' @return a list of a dataset (list of a matrix mutation count and a matrix of read depth), a vector z of group assignment and a list of parameters
simu_hidden_RW <- function (n, time, date.origin = NULL, param, lambda) {
 
  pi = param$pi
  mu = param$mu
  s = param$s
  alpha = param$alpha
  beta = param$beta
  sd0 = param$sd0
  sd = param$sd
  
  # number of time points
  nT = length(time)
  delta = diff(time) # difference between to subsequent time points
  
  # number of non-neutral groups 
  K = length(pi) - 1
  
  # group assignment
  z = NULL; while(length(unique(z))!=K+1) {z = sample(0:K, n, replace = TRUE, prob = pi)}
  
  # Randam walk
  y = matrix(NA, n, nT)
  # Y0 conditional on Z=0
  u = rbeta(sum(z==0), alpha, beta)
  y[z==0, 1] = log(u/(1-u)) 
  # Y0 conditional on Z=k, k\neq 0
  for (k in 1:K) {
    y[z==k, 1] = rnorm(sum(z==k), mu[k], sd0)
  }
  # fill for t \neq 0
  for (t in 2:nT) {
    y[,t] = rnorm(n, y[,t-1] + c(0.00, s)[z+1]*delta[t-1], sd*sqrt(delta[t-1]))
  }
  
  # read depth
  d = matrix(rpois(n * nT, lambda = lambda), n, nT)
  # mutation count
  x = matrix(rbinom(n*nT, d, exp(y)/(1+exp(y))), nrow = n, ncol = nT, byrow= FALSE)
  
  # Convert colnames x and d to be in line with real data
  colnames(x) = colnames(d) = as.character(as.Date (ifelse(is.null(date.origin), "2019-01-01", date.origin)) + time)
  data = list(x = x, d = d)
  
  return(list(data = list(x=x, d=d), z = z))
}



