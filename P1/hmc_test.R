K <- function(p){
    return(sum(t(p) %*% p) / 2)
}

# Auxilliary variables for visualisation
x=seq(-5,5,by=0.1)
y=seq(-5,5,by=0.1)

z1 = matrix(0, ncol=101, nrow=101)
z2 = matrix(0, ncol=101, nrow=101)
z3 = matrix(0, ncol=101, nrow=101)


##############################################
# 1. Standard Gaussian with Correlation

# Parameters
Sigma = matrix(c(1,0.9,0.9,1), nrow=2)
SigmaInv = solve(Sigma)
SigmaSqrt = sqrt(det(Sigma))

# Density
myGaussian <- function(x){
    f = 1/(2*pi*SigmaSqrt)*exp(-0.5*t(x)%*%SigmaInv%*%x)
    return(f)
}


##############################################
# 2. Multimodal

# Parameters
mu1 = c(-1.5,-1.5)
Sigma1 = matrix(c(1,0,0,1), ncol=2)
Sigma1Inv = solve(Sigma1)
Sigma1Sqrt = sqrt(det(Sigma1))

mu2 = c(1.5,1.5)
Sigma2 = matrix(c(1,0,0,1), ncol=2)
Sigma2Inv = solve(Sigma2)
Sigma2Sqrt = sqrt(det(Sigma2))

mu3 = c(-2,2)
Sigma3 = matrix(c(0.8,0,0,0.8), ncol=2)
Sigma3Inv = solve(Sigma3)
Sigma3Sqrt = sqrt(det(Sigma3))

# Density
myMultimodal <- function(x){
    f = 1/3 * 1/(2*pi*Sigma1Sqrt)*exp(-0.5*t(x-mu1)%*%Sigma1Inv%*%(x-mu1))
    f = f + 1/3 * 1/(2*pi*Sigma2Sqrt)*exp(-0.5*t(x-mu2)%*%Sigma2Inv%*%(x-mu2))
    f = f + 1/3 * 1/(2*pi*Sigma3Sqrt)*exp(-0.5*t(x-mu3)%*%Sigma3Inv%*%(x-mu3))
    return(f)
}

##############################################
# 3. Volcano

# Density
myVolcano <- function(x){
    f = 1/(2*pi)*exp(-0.5*t(x)%*%x)*(t(x)%*%x+0.25)
    return(f)
}


mcmc <- function(x0, n, proposal_func, acceptance_func, dens) {
    
    # object to save MCMC samples
    x <- matrix(0, nrow=n, ncol=2)
    alphas <- rep(0,n)
    # initialisation
    x_old <- x0
    # generate a vector of n uniform distributed random variables
    u <- runif(n)
    # go through n itearations
    for (i in 1:n) {
        # get a new proposal for x
        x_prop <- proposal_func(x_old)
        # compute the acceptance prob alpha
        alpha <- acceptance_func(x_prop, x_old, dens)
        alphas[i] = alpha
        # decide wether to accept or reject
        if (u[i] < alpha) {
            x_new <- x_prop
        } else {
            x_new <- x_old
        }
        # save the new sample
        x[i,] <- x_new
        x_old <- x_new
    }
    return(list("x"=x,"alphas"=alphas))
}

HMC_acceptance <- function(x_prop, x_old, dens){
    U_old <- -log(dens(x_old))
    U_prop <- -log(dens(x_prop))
    alpha = min(1, exp(U_old - U_prop + K_old - K_prop))
    return(alpha)
}

HMC_proposal <- function(x_old){
    x <- x_old
    p <- rnorm(length(x), 0, 1)
    p_old <- p
    res_leapfrog = leapfrog(p, x, delta, L, dU)
    p = - res_leapfrog[, 1]
    x = res_leapfrog[, 2]
    K_old = K(p_old)
    K_prop = K(p)
    # print(K_old)
    # print(K_prop)
    return(x)
}

leapfrog <- function(p, x, delta, L, dU){
    p <- p - delta * dU(x) / 2
    for (i in 1:L){
        x <- x + delta * p
        if (i != L) p <- p - delta * dU(x)
    }
    p <- p - delta * dU(x) / 2
    # print(p)
    # print(x)
    return(cbind(p, x))
}

# leapfrog <- function(p, x, delta, L, dU_test){
#     p <- p - delta * dU_test(x) / 2
#     for (i in 1:L){
#         x <- x + delta * p
#         if (i != L) p <- p - delta * dU_test(x)
#     }
#     p <- p - delta * dU_test(x) / 2
#     print(p)
#     print(x)
#     return(cbind(p, x))
# }

# dU <- function(x_old){
#     du = rep(0,length(x_old))
#     h = 1e-3
#     e1 = c(1,0)
#     e2 = c(0,1)
#     du[1] = (-log(dens(x_old-h*e1)) + log(dens(x_old+h*e1)))/(2*h)
#     du[2] = (-log(dens(x_old-h*e2)) + log(dens(x_old+h*e2)))/(2*h)
#     return(du)
# }
U <- function(x){
    return(-log(dens(x)))
}

dU <- function(x_old){
    du = grad(U, x_old)
    if (any(is.nan(du))){
        du = c(0, 0)
    }
    return(du)
}


x0 = c(0,0)
n = 10000
delta = 0.3
L = 20
# dens = myGaussian
# dens = myMultimodal
dens = myVolcano
K_old = 0
K_prop = 0

trace0 = mcmc(x0, n, HMC_proposal, HMC_acceptance, dens)
MCMCplots(trace0)


# test = leapfrog(t(t(c(2, 1))), t(t(c(2, 1))), 1, 10, sin)

MCMCplots <- function(trace, surface=FALSE){
    # Generating empirical densitty
    kd = MASS::kde2d(x=trace$x[,1], y=trace$x[,2], n=101, lims=c(-5,5,-5,5))
    
    # Nice 3D empirical density plots 
    # WARNING: VERY SLOW
    if (surface==TRUE){
        library(plotly)
        plot_ly(x=kd$x, y=kd$y, z=kd$z, type="surface")
    }
    
    # 2D empirical density plots
    image( x=kd$x, y=kd$y, z=kd$z,
           asp=1, main="Density",
           xlab="x", ylab="y", xlim=c(-5,5), ylim=c(-5,5))
    
    # Trace plots
    # plot(trace$x[,1], type="l", ylab="x1", 
    #      main=paste("Trace of x1 (tuning parameter",sigma,")"))
    # plot(trace$x[,2], type="l", ylab="x1", 
    #      main=paste("Trace of x2 (tuning parameter",sigma,")"))
    
    # Autocorrelation
    acf(trace$x[,1], main="Trace plot for x1")
    acf(trace$x[,2], main="Trace plot for x2")
    
    # Mean acceptance rate
    print(paste("The mean acceptance rate is", sum(trace$alphas)/length(trace$alphas)))
}




# setup
# delta = 0.3
# nSamples = 10000
# L = 20
# 
# # initialisation
# x = matrix(0, 2, nSamples)
# x[, 1] = c(0, 6)
# t = 1
# while(t < nSamples){
#     t = t + 1
#     x[, t] = HMC(U = U1, dU = dU1, delta = delta, L = L, current_q = t(t(x[, t-1])))
#     # x[, t] = HMC(U = U2, dU = dU2, delta = delta, L = L, current_q = t(t(x[, t-1])))
#     # x[, t] = HMC(U = U2, dU = dU2, delta = delta, L = L, current_q = t(t(x[, t-1])))
# }
# 
# plot(x[1, ])
# plot(x[1, ], x[2, ])
# par(new = TRUE)
# plot(x[1, 1:500], x[2, 1:500], type = "b", pch = 19, col = 'red')
# library(ggplot2)
# library(gridExtra)
# # plot(x[1, 1:1000], type = "l")
# hist_top <- ggplot() + geom_histogram(aes(x[1, ]))
# empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + 
#          theme(axis.ticks = element_blank(), 
#                panel.background = element_blank(), 
#                axis.text.x = element_blank(), axis.text.y = element_blank(), 
#                axis.title.x = element_blank(), axis.title.y = element_blank())
# scatter <- ggplot() + geom_point(aes(x[1, ], x[2, ]))
# hist_right <- ggplot() + geom_histogram(aes(x[2, ])) + coord_flip()
# 
# grid.arrange(hist_top, empty, scatter, hist_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
# 
# 
# 
# 
# 
