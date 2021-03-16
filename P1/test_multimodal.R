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

z2 = matrix(0, ncol=101, nrow=101)
x=seq(-5,5,by=0.1)
y=seq(-5,5,by=0.1)
# Fill z
for (i in 1:101){
    for (j in 1:101){
        z2[i,j] = myMultimodal(c(x[i],y[j]))
    }
}

# Visualise
image( x=x,
       y=y,
       z=z2,
       asp=1,
       main="Density of Multimodal",
       xlab="x", ylab="y")


# ====================== # ======================# ======================# ======================
# back up codes

# # Gaussian
# U1 <- function(x){
#     U = t(x) %*% solve(matrix(c(1,0.9,0.9,1), nrow=2)) %*% x
#     return(U)
# }
# 
# dU1 <- function(x){
#     dU = t(x) %*% solve(matrix(c(1,0.9,0.9,1), nrow=2))
#     return(dU)
# }

# # Multimodal
# U2 <- function(x){
#     mu1 = c(-1.5,-1.5)
#     mu2 = c(1.5,1.5)
#     mu3 = c(-2,2)
#     U = t(x - mu1) %*% solve(matrix(c(1,0,0,1), ncol=2)) %*% (x - mu1) + 
#         t(x - mu2) %*% solve(matrix(c(1,0,0,1), ncol=2)) %*% (x - mu2) + 
#         t(x - mu3) %*% solve(matrix(c(0.8,0,0,0.8), ncol=2)) %*% (x - mu3)
#     return(U)
# }
# 
# dU2 <- function(x){
#     mu1 = c(-1.5,-1.5)
#     mu2 = c(1.5,1.5)
#     mu3 = c(-2,2)
#     dU = t(x - mu1) %*% solve(matrix(c(1,0,0,1), nrow=2)) + 
#          t(x - mu2) %*% solve(matrix(c(1,0,0,1), nrow=2)) + 
#          t(x - mu3) %*% solve(matrix(c(0.8,0,0,0.8), nrow=2)) 
#     return(dU)
# }

# # Volcano
# U3 <- function(x){
#     U = -1/2 * t(x) %*% x * log(t(x) %*% x + 0.25)
#     return(U)
# }
# dU3 <- function(x){
#     dU = -t(x) %*% log(x %*% t(x) + 0.25) - 1/2 * t(x) %*% x %*% t(x) %*% solve(x %*% t(x) + 0.25)
#     return(dU)
# }

U1 <- function(x){
    return(-log(myGaussian(x)))
}

dU1 <- function(x){
    return(dU(x, myGaussian))
}

U2 <- function(x){
    return(-log(myMultimodal(x)))
}

dU2 <- function(x){
    return(dU(x, myMultimodal))
}

U3 <- function(x){
    return(-log(myVolcano(x)))
}

dU3 <- function(x){
    return(dU(x, myVolcano))
}




# U is a function that returns the potential energy given q
# grad_U returns the respective partial derivatives
# epsilon stepsize
# L number of leapfrog steps
# current_q current position

# kinetic energy is assumed to be sum(p^2/2) (mass == 1)
HMC <- function (U, dU, delta, L, current_q) {
    q <- current_q # x <- x_old
    # independent standard normal variates
    p <- t(rnorm(length(q), 0, 1))
    # Make a half step for momentum at the beginning
    current_p <- p
    # Alternate full steps for position and momentum
    p <- p - delta * dU(q) / 2
    for (i in 1:L) {
        # Make a full step for the position
        q <- q + delta * t(p)
        # Make a full step for the momentum, except at end of trajectory
        if (i != L) p <- p - delta * dU(q)
    }
    # Make a half step for momentum at the end
    p <- p - delta * dU(q) / 2
    # Negate momentum at end of trajectory to make the proposal symmetric
    p <- -p
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U <- U(current_q)
    current_K <- K(current_p)  #sum(current_p^2) / 2
    proposed_U <- U(q)
    proposed_K <- K(p) #sum(p^2) / 2
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
        return (q)  # accept
    } else {
        return (current_q)  # reject
    }
}





