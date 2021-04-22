# 
# Matern_cov <- function(sigma, phi, t){
#     # param sigma: scaling coef
#     # param eta: range coef
#     # param t: distance matrix
#     # return: matern covariance
#     return(sigma ^ 2 * (1 + phi * t) * exp(-phi * t))
# }
# 
# plotf2d <- function(v, string){
#   vv <- v
#   dim(vv) <- c(n1, n2)
#   levelplot(vv, col.regions = coul, main = string)
# }
# 
# # Setup the grid
# n1 = 25 # number of grid points along east direction
# n2 = 25 # number of grid points along north direction
# n = n1 * n2 # total number of grid points
# 
# dn1 = 1/n1
# dn2 = 1/n2
# sites1 = array(seq(0, 1, dn1), c(n1, 1))
# sites2 = array(seq(0, 1, dn2), c(n2, 1))
# 
# ww1 = rep(1, n1)
# ww2 = rep(1, n2)
# sites1m = sites1 %*% t(ww1) # sites1m is the matrix version of sites1
# sites2m = ww2 %*% t(sites2)
# 
# sites1v = matrix(sites1m, nrow = n, ncol = 1)
# sites2v = matrix(sites2m, nrow = n, ncol = 1)
# 
# # plot(sites1v, sites2v)
# 
# # Compute the distance matrix
# ddE = sites1v %*% matrix(rep(1, n), nrow = 1, ncol = n) - matrix(rep(1, n), nrow = n, ncol = 1) %*% t(sites1v)
# dd2E = ddE * ddE
# ddN = sites2v %*% matrix(rep(1, n), nrow = 1, ncol = n) - matrix(rep(1, n), nrow = n, ncol = 1) %*% t(sites2v)
# dd2N = ddN * ddN
# t = sqrt(dd2E + dd2N)
# levelplot(t, col.regions = coul, main = "Distance matrix")
# 
# # Simulate the initial random field
# alpha = 1.0 # beta as in regression model
# sigma = 1.0  # scaling coef in matern kernel
# phi = 10 # range coef in matern kernel
# # eta = 10 # range coef in matern kernel
# tau = .05 # iid noise
# 
# beta1 = -alpha
# beta2 = alpha
# beta3 = alpha
# 
# BETA_TRUE = rbind(beta1, beta2, beta3)
# THETA_TRUE = rbind(sigma, phi, tau)
# 
# Sigma = Matern_cov(sigma, phi, t)  # matern covariance
# 
# L = t(chol(Sigma)) # lower L
# x = L %*% rnorm(n) # sample from zero mean random variables
# 
# H = array(c(rep(1, n), sites1v, sites2v), dim = c(n, 3)) # design matrix
# mu_prior = H %*% BETA_TRUE
# plotf2d(mu_prior, "Prior mean")
# 
# mu_real = mu_prior + x
# plotf2d(mu_real, "Realisation of the fields")
# 
# # sampling from realisations in the regular grids
# M = 200
# Fmatrix = matrix(0, M, n)
# ind = sample(n, size = M, replace = FALSE)
# for (i in c(1:M)){
#   Fmatrix[i, ind[i]] = TRUE
# }
# G = Fmatrix %*% H
# y_sampled = Fmatrix %*% mu_real + tau * rnorm(M, 1)
# x_ind = sites1v[ind]
# y_ind = sites2v[ind]
# 
# plot(x_ind, y_ind, cex = abs(y_sampled), main = "Random samples in the grid, circle size indicates the relative value")
# 
# 
# # sampling randomly without grids
# while (TRUE){
#   x_ind = runif(M)
#   y_ind = runif(M)
#   if ((unique(x_ind) == x_ind) && (unique(y_ind) == y_ind)){
#     dim(x_ind) <- c(M, 1)
#     dim(y_ind) <- c(M, 1)
#     break;
#   }
# }
# 
# ddE_sample = x_ind %*% matrix(rep(1, M), nrow = 1, ncol = M) - matrix(rep(1, M), nrow = M, ncol = 1) %*% t(x_ind)
# dd2E_sample = ddE_sample * ddE_sample
# ddN_sample = y_ind %*% matrix(rep(1, M), nrow = 1, ncol = M) - matrix(rep(1, M), nrow = M, ncol = 1) %*% t(y_ind)
# dd2N_sample = ddN_sample * ddN_sample
# tsample = sqrt(dd2E_sample + dd2N_sample)
# Sigma_sample = Matern_cov(sigma, phi, tsample)
# Lsample = t(chol(Sigma_sample))
# xsample = Lsample %*% rnorm(M)
# Hsample = array(c(rep(1, M), x_ind, y_ind), dim = c(M, 3)) # sampling design matrix
# mu_prior_sample = Hsample %*% BETA_TRUE
# mu_real_sample = mu_prior_sample + xsample

# plot(x_ind, y_ind, cex = abs(mu_real_sample), main = "Random samples in the square, circle size indicates the relative value")
# y_sampled = mu_real_sample + tau * rnorm(M, 1)




### 2.1.1 Simulation using a grid

To simulate from the grid discretisation, the $\boldsymbol{F}$ matrix is defined, which contains the information on where to sample. \ref{fig:SamplingMatrix} is a 25 $\times$ 25 demonstrates the sampling matrix for one step, i.e., the red indicates where it has sampled, therefore, it becomes 1 in the matrix, the rest are zeros. By having multiple steps, and vectorise each matrix to a long vector which is a m $\times$ N matrix, N is 625. The red grids are regions where it have sampled as shown in \ref{fig:FMatrix}. They are only illustration for the purpose of understanding of how this grid approach works. It has the benefit of improving the numerical stabilities since one grid can only contain one value. Otherwise, very close value might induce singularity in later inverse computations. 

```{r Sampling matrix, echo=FALSE, fig.cap="\\label{fig:SamplingMatrix}Sampling design matrix in one step", out.width = '100%'}
knitr::include_graphics("F.pdf")
```
```{r F matrix, echo=FALSE, fig.cap="\\label{fig:FMatrix}F matrix after 10 steps", out.width = '100%'}
knitr::include_graphics("Fv.pdf")
```